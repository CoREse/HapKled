#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <vector>
#include <string>
#include <mutex>
#include <unistd.h>
//#include <omp.h>
#include "htslib/htslib/bgzf.h"
#include <thread>
#include <stdexcept>
#include "input.h"
#include "signature.h"
#include <filesystem>
#include "contig.h"
#include "crelib/crelib.h"
#include "StatsManager.hpp"
#include <math.h>
#include <algorithm>
#include "htslib/thread_pool.h"

using namespace std;
using namespace cre;

const int ReadThreadN=1;//old mechanism that leaves this number of threads handling reads. not good. obsoleted. use htslib's thread mechanisms instead

float MedianInsertionSize=481.2;
float ISSD=36.4;
//bool DynamicDRP=false;
//int DynamicCount=1000;
int RDWindowSize=100;

mutex TopLock, WriteLock, StdinLock;

struct HandleBrMutex
{
	pthread_mutex_t m_AllPrimarySeg;
	pthread_mutex_t m_Cov;
	pthread_mutex_t m_Sig[NumberOfSVTypes];
	HandleBrMutex()
	{
		pthread_mutex_init(&m_AllPrimarySeg,NULL);
		pthread_mutex_init(&m_Cov,NULL);
		for (int i=0;i<NumberOfSVTypes;++i)
		{
			pthread_mutex_init(m_Sig+i,NULL);
		}
	}
};

int getTechFromReads(bam_hdr_t *Header, htsFile* SamFile)
{
	int CheckN=10;
	bam1_t *br=bam_init1();
	for (int i=0;i<CheckN;++i)
	{
		sam_read1(SamFile,Header,br);
		if (bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br))>1000) return 0;
	}
	return 1;
}

int getTech(const char * ReferenceFileName, const char * BamFileName)//return 0: SMRT, 1:NGS, -1: others
{
	htsFile* SamFile;//the file of BAM/CRAM
	bam_hdr_t* Header;//header for BAM/CRAM file
	hts_idx_t* BamIndex;
	vector<Stats> AllStats;
	const char * SampleFileName=BamFileName;
	SamFile = hts_open(SampleFileName, "rb");
	BamIndex= sam_index_load(SamFile,SampleFileName);
	if (BamIndex==NULL)
	{
		hts_close(SamFile);
		die("Please index the samfile(%s) first!",SampleFileName);
	}
	//set reference file
	if (NULL != ReferenceFileName)
	{
		char referenceFilenameIndex[128];
		strcpy(referenceFilenameIndex, ReferenceFileName);
		strcat(referenceFilenameIndex, ".fai");
		int ret = hts_set_fai_filename(SamFile, referenceFilenameIndex);
	}
	Header = sam_hdr_read(SamFile);

	int Tech;//0: SMRT, 1:NGS
	kstring_t ks=KS_INITIALIZE;
	ks_resize(&ks,1000);
	char PL[1000];
	sam_hdr_find_line_pos(Header,"RG",0,&ks);
	int splitn;
	int * splitoffsets=ksplit(&ks,'\t',&splitn);
	for (int i=0;i<splitn;++i)
	{
		if (memcmp(ks.s+splitoffsets[i],"PL",2)==0)
		{
			int End=ks_len(&ks);
			if (i<splitn-1) End=splitoffsets[i+1];
			memcpy(PL,ks.s+splitoffsets[i]+3,End-splitoffsets[i]-3);
			PL[End-splitoffsets[i]-3]='\0';
		}
	}
	ks_free(&ks);
	if (strcmp(PL,"ILLUMINA")==0)
	{
		Tech=1;
	}
	else if (strcmp(PL,"PACBIO")==0 || strcmp(PL,"ONT")==0)
	{
		Tech=0;
	}
	else
	{
		Tech=getTechFromReads(Header,SamFile);
	}
	hts_close(SamFile);
	return Tech;
}

vector<int> getAllTechs(Arguments & Args)
{
	vector<int> AllTechs;
	for (int i=0;i<Args.BamFileNames.size();++i)
	{
		const char * SampleFileName=Args.BamFileNames[i];
		AllTechs.push_back(getTech(Args.ReferenceFileName,Args.BamFileNames[i]));
	}
	return AllTechs;
}

//get soft and hard clip length sum
float getClipLength(bam1_t * br)
{
	uint32_t* cigars=bam_get_cigar(br);
	float ClipLength=0;
	for (int i=0;i<br->core.n_cigar;++i)
	{
		if ((cigars[i]&0xf)==4 || (cigars[i]&0xf)==5) ClipLength+=cigars[i]>>4;
	}
	return ClipLength;
}

//get query length without head and tail clips
#define brGetClippedQlen(br) (getClippedQLen(br->core.n_cigar,bam_get_cigar(br)))
//get query length without head and tail clips
int getClippedQLen(uint32_t CIGARN, uint32_t* CIGARD)
{
	int Start=0,N=CIGARN;
	if (bam_cigar_op(CIGARD[0])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[0])==BAM_CHARD_CLIP)
	{
		Start=1;
		--N;
	}
	if (bam_cigar_op(CIGARD[N])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[N])==BAM_CHARD_CLIP)
	{
		--N;
	}
	return bam_cigar2qlen(N,CIGARD+Start);
}

//get read length(with all clipped length)
#define brGetReadLength(br) (getReadLength(br->core.n_cigar,bam_get_cigar(br)))
//get read length(with all clipped length)
int getReadLength(uint32_t CIGARN, uint32_t* CIGARD)//included hard clip
{
	int Start=0,N=CIGARN;
	int Clipped=0;
	if (bam_cigar_op(CIGARD[0])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[0])==BAM_CHARD_CLIP)
	{
		Start=1;
		--N;
		Clipped+=bam_cigar_oplen(CIGARD[0]);
	}
	if (bam_cigar_op(CIGARD[CIGARN-1])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[CIGARN-1])==BAM_CHARD_CLIP)
	{
		--N;
		Clipped+=bam_cigar_oplen(CIGARD[CIGARN-1]);
	}
	return bam_cigar2qlen(N,CIGARD+Start)+Clipped;
}

void swap(int &a, int &b)
{
	int c=a;
	a=b;
	b=c;
}

struct Alignment
{
	int Cid;//Contig id of mapping position
	int Pos;
	int Length;//Length in Reference
	int End;
	int Strand;//0: -, 1: +
	int MapQ;
	int NM;
	vector<uint32_t> CIGAR;
	int InnerPos;//for same orient segments, just use innerpos, forward pos if for sorting of mixed orientations
	int InnerLength;
	int InnerEnd;
	int ForwardPos;//for inner sort
	int ForwardEnd;
	Alignment(int Cid, int Pos, int Strand, uint32_t* CIGARD, uint32_t CIGARN, int mapQ, int NM)
	:Cid(Cid), Pos(Pos), Strand(Strand), MapQ(mapQ), NM(NM)
	{
		construct(Pos,Strand,CIGARD,CIGARN);
	}
	void construct(int Pos, int Strand, uint32_t* CIGARD, uint32_t CIGARN)
	{
		this->Pos=Pos;
		this->Strand=Strand;
		Length=bam_cigar2rlen(CIGARN, CIGARD);
		End=Pos+Length;
		if (bam_cigar_opchr(*CIGARD)=='H' || bam_cigar_opchr(*CIGARD)=='S') InnerPos=bam_cigar_oplen(*CIGARD);
		else InnerPos=0;
		InnerLength=getClippedQLen(CIGARN,CIGARD);
		InnerEnd=InnerPos+InnerLength;
		if (Strand==0)
		{
			int ReadLength=getReadLength(CIGARN,CIGARD);
			ForwardPos=ReadLength-InnerEnd;
			ForwardEnd=ReadLength-InnerPos;
			InnerPos=ForwardPos;
			InnerEnd=ForwardEnd;
		}
		else
		{
			ForwardPos=InnerPos;
			ForwardEnd=InnerEnd;
		}
		for (int i =0;i<CIGARN;++i) CIGAR.push_back(CIGARD[i]);
	}
	Alignment(int cid, int Pos, int Strand, const char * CIGARS, int MapQ, int NM)
	{
		this->Cid=cid;
		this->MapQ=MapQ;
		this->NM=NM;
		size_t CIGARN=size_t(strlen(CIGARS));//sam_parse_cigar will reallocate if not enough
		uint32_t * ca=(uint32_t*) malloc(CIGARN*sizeof(uint32_t));
		CIGARN=sam_parse_cigar(CIGARS,NULL,&ca,&CIGARN);
		construct(Pos,Strand,ca,CIGARN);
		free(ca);
	}
	bool operator<(const Alignment &other) const
	{
		return this->ForwardPos<other.ForwardPos;
		// return this->InnerPos<other.InnerPos;
	}
	double getQuality()
	{
		return 1.0-((double)(NM))/((double)InnerLength);//InnerLength-NM to estimate AS
	}
};

string cigar2string(uint32_t* CIGARD, uint32_t CIGARN)
{
	string result="";
	for (int i=0;i<CIGARN;++i)
	{
		result+=to_string(bam_cigar_oplen(CIGARD[i]));
		result+=bam_cigar_opchr(CIGARD[i]);
	}
	return result;
}

void searchDelFromAligns(bam1_t *br, Contig& TheContig, vector<Alignment> &Aligns, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, HandleBrMutex *Mut, Arguments & Args)
{
	vector<Signature> &Signatures=TypeSignatures[TheContig.ID][0];
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=i+1;j<Aligns.size();++j)
		{
			if (Aligns[i].Strand!=Aligns[j].Strand || Aligns[i].Cid != Aligns[j].Cid)
			{
				continue;
				// if (i+1<Aligns.size())
				// {
				// 	if (Aligns[i-1].Strand!=Aligns[i+1].Strand) continue;
				// 	i2=i+1;
				// }
				// else continue;
			}
			if (abs(Aligns[i].Pos-Aligns[j].Pos)>1000000) continue;
			int FormerI=i, LatterI=j;
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			// if (Aligns[i-1].End<Aligns[i].Pos && Aligns[i].Pos-Aligns[i-1].End-(Aligns[i].InnerPos-Aligns[i-1].InnerEnd)>=Args.MinSVLen) Signatures.push_back(Signature(2,Tech,0,Aligns[i-1].End,Aligns[i].Pos,bam_get_qname(br),br->core.qual));
			if (Aligns[FormerI].End<Aligns[LatterI].Pos && Aligns[LatterI].Pos-Aligns[FormerI].End-(Aligns[LatterI].InnerPos-Aligns[FormerI].InnerEnd)>=Args.MinSVLen)
			{
				// int End=Aligns[FormerI].End+Aligns[LatterI].Pos-Aligns[FormerI].End-(Aligns[LatterI].InnerPos-Aligns[FormerI].InnerEnd);
				int End=Aligns[LatterI].Pos;
				double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
				Signature TempSignature(2,Tech,0,Aligns[FormerI].End,End,bam_get_qname(br),Quality);
				pthread_mutex_lock(&Mut->m_Sig[0]);
				TypeSignatures[Aligns[i].Cid][0].push_back(TempSignature);
				pthread_mutex_unlock(&Mut->m_Sig[0]);
			}
			break;//Get the first one only. Prevent duplicate
		}
	}
}

void searchInsFromAligns(bam1_t *br,Contig& TheContig,vector<Alignment> &Aligns, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, HandleBrMutex *Mut, Arguments & Args)
{
	// vector<Signature> &Signatures=TypeSignatures[TheContig.ID][1];
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=i+1;j<Aligns.size();++j)
		{
			if (Aligns[i].Strand!=Aligns[j].Strand) continue;
			if (abs(Aligns[i].Pos-Aligns[j].Pos)>1000000) continue;
			int FormerI=i, LatterI=j;
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			int Gap=(Aligns[LatterI].InnerPos-Aligns[FormerI].InnerEnd)-(Aligns[LatterI].Pos+10-Aligns[FormerI].End);
			if (Aligns[FormerI].End<Aligns[LatterI].Pos+Args.InsClipTolerance && Gap>=Args.MinSVLen && Gap<Args.InsMaxGapSize)
			{
				double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
				Signature TempSignature(2,Tech,1,(Aligns[FormerI].End+Aligns[LatterI].Pos)/2,MIN(TheContig.Size-1,(Aligns[FormerI].End+Aligns[LatterI].Pos)/2+Gap),bam_get_qname(br),Quality);
				pthread_mutex_lock(&Mut->m_Sig[1]);
				TypeSignatures[Aligns[FormerI].Cid][1].push_back(TempSignature);
				// Signatures.push_back(TempSignature);
				pthread_mutex_unlock(&Mut->m_Sig[1]);
			}
			break;
		}
	}
}

void searchDupFromAligns(bam1_t *br,Contig& TheContig,vector<Alignment> &Aligns, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, HandleBrMutex *Mut, Arguments & Args)
{
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=0;j<Aligns.size();++j)
		{
			if (Aligns[i].Strand!=Aligns[j].Strand) continue;
			if (abs(Aligns[i].Pos-Aligns[j].Pos)>1000000) continue;
			// vector<Segment> v;
			int FormerI=i, LatterI=j;
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			if (Aligns[LatterI].Pos<Aligns[FormerI].End)
			{
				int Dup=0;
				if (Aligns[FormerI].End-Aligns[LatterI].Pos>Aligns[LatterI].Length)
				{
					Dup=1;
				}
				else if (Aligns[FormerI].End-Aligns[LatterI].Pos>Aligns[FormerI].Length)
				{
					Dup=1;
				}
				else
				{
					if (Aligns[FormerI].End-Aligns[LatterI].Pos>=Args.MinSVLen)
					{
						Dup=1;
					}
				}
				if (Dup)
				{
					double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
					// v.push_back(Segment(Aligns[i-1].Pos,Aligns[i-1].End));
					// v.push_back(Segment(Aligns[i].Pos,Aligns[i].End));
					Signature TempSignature(2,Tech,2,Aligns[LatterI].Pos,Aligns[FormerI].End,bam_get_qname(br),Quality);
					pthread_mutex_lock(&Mut->m_Sig[2]);
					TypeSignatures[Aligns[FormerI].Cid][2].push_back(TempSignature);
					pthread_mutex_unlock(&Mut->m_Sig[2]);
				}
			}
			//if (Aligns[i-1].End<Aligns[i].Pos && Aligns[i].Pos-Aligns[i-1].End-(Aligns[i].InnerPos-Aligns[i-1].InnerEnd)>=50) Signatures.push_back(Signature(2,Tech,2,Aligns[i-1].End,Aligns[i].Pos,bam_get_qname(br),br->core.qual));
		}
		break;
	}
}

bool continuous(const Alignment& Former, const Alignment& Latter, unsigned Endurance)
{
	// return true;
	// if(abs(Latter.Pos-Former.End)<Endurance && abs(Latter.InnerPos-Former.InnerEnd)<Endurance) return true;
	if (abs(Latter.Pos-Former.Pos)>1000000) return false;
	if (Latter.Pos+Endurance>=Former.End && Latter.Pos>Former.Pos) return true;
	return false;
}

void searchInvFromAligns(bam1_t *br,Contig& TheContig,vector<Alignment> &Aligns, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, HandleBrMutex *Mut, Arguments & Args)
{
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=i+1;j<=i+1;++j)
		{
			if (Aligns[i].Strand==Aligns[j].Strand) continue;
			if (Aligns[i].Cid!=Aligns[j].Cid) continue;
			int InvClipEndurance=1000;
			int FormerI=i, LatterI=j;//k is always the middle one, former and latter is for ref position
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			if (continuous(Aligns[FormerI],Aligns[LatterI],InvClipEndurance))
			{
				int k=0;
				for (k=j+1;k<Aligns.size();++k)
				{
					if (Aligns[i].Cid!=Aligns[k].Cid) continue;
					if (Aligns[i].Strand==0)
					{
						if (Aligns[j].Strand!=Aligns[k].Strand && continuous(Aligns[k],Aligns[j],InvClipEndurance)) break;
					}
					else
					{
						if (Aligns[j].Strand!=Aligns[k].Strand && continuous(Aligns[j],Aligns[k],InvClipEndurance)) break;
					}
				}
				if (k!=Aligns.size())
				{
					double Quality=Aligns[i].getQuality()+Aligns[j].getQuality()+Aligns[k].getQuality();
					Quality/=3.0;
					Signature TempSignature(2,Tech,3,Aligns[j].Pos,Aligns[j].End,bam_get_qname(br),Quality);
					TempSignature.setInvLeft(true);
					TempSignature.setInvRight(true);
					pthread_mutex_lock(&Mut->m_Sig[3]);
					TypeSignatures[Aligns[i].Cid][3].push_back(TempSignature);
					pthread_mutex_unlock(&Mut->m_Sig[3]);
					++i;
				}
				else
				{
					double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
					Signature Temp1Signature(2,Tech,3,Aligns[FormerI].Pos,Aligns[FormerI].End,bam_get_qname(br),Quality);
					Temp1Signature.setInvRight(true);
					Signature Temp2Signature(2,Tech,3,Aligns[LatterI].Pos,Aligns[LatterI].End,bam_get_qname(br),Quality);
					Temp2Signature.setInvLeft(true);
					pthread_mutex_lock(&Mut->m_Sig[3]);
					TypeSignatures[Aligns[FormerI].Cid][3].push_back(Temp1Signature);
					TypeSignatures[Aligns[FormerI].Cid][3].push_back(Temp2Signature);
					pthread_mutex_unlock(&Mut->m_Sig[3]);
				}
				break;
			}
		}
	}
}

inline void statCoverage(int Begin, int End, double *CoverageWindows, Contig & TheContig, Arguments &Args)
{
	if (End>=TheContig.Size) End=TheContig.Size-1;
	if (End<=Begin) return;
	int WBegin=Begin/Args.CoverageWindowSize;
	int WEnd=End/Args.CoverageWindowSize+1;
	for (int i=WBegin+1;i<WEnd-1;++i) CoverageWindows[i]+=1;
	double FirstPortion=((double)(((WBegin+1)*Args.CoverageWindowSize)-Begin))/((double)Args.CoverageWindowSize);
	CoverageWindows[WBegin]+=FirstPortion;
	if (WEnd>WBegin+1)
	{
		double LastPortion=((double)(End+1-(WEnd-1)*Args.CoverageWindowSize))/((double)Args.CoverageWindowSize);
		CoverageWindows[WEnd-1]+=LastPortion;
	}
}
inline void statCoverageCigar(bam1_t * br, double *CoverageWindows, Contig & TheContig, HandleBrMutex *Mut, Arguments &Args)
{
	// printf("%s %d %d\n", bam_get_qname(br),br->core.pos, br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br)));
	if (br->core.qual<Args.MinMappingQuality) return;
	int Begin=br->core.pos;
	int End=Begin+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br));
	if (End>=TheContig.Size) End=TheContig.Size-1;
	if (End<=Begin) return;
	int WBegin=Begin/Args.CoverageWindowSize;
	int WEnd=End/Args.CoverageWindowSize+1;
	double FirstPortion=((double)(((WBegin+1)*Args.CoverageWindowSize)-Begin))/((double)Args.CoverageWindowSize);
	pthread_mutex_lock(&Mut->m_Cov);
	for (int i=WBegin+1;i<WEnd-1;++i) CoverageWindows[i]+=1;
	CoverageWindows[WBegin]+=FirstPortion;
	pthread_mutex_unlock(&Mut->m_Cov);
	if (WEnd>WBegin+1)
	{
		double LastPortion=((double)(End+1-(WEnd-1)*Args.CoverageWindowSize))/((double)Args.CoverageWindowSize);
		pthread_mutex_lock(&Mut->m_Cov);
		CoverageWindows[WEnd-1]+=LastPortion;
		pthread_mutex_unlock(&Mut->m_Cov);
	}
}

void dealClipConflicts(vector<Alignment> &Aligns, Arguments & Args)
{
	int Endurance=100;
	for (int i=1;i<Aligns.size();++i)
	{
		if (Aligns[i].Strand!=Aligns[i-1].Strand) continue;
		int FormerI=i-1, LatterI=i;
		if (Aligns[i-1].Strand==0)
		{
			FormerI=i;
			LatterI=i-1;
		}
		if (Aligns[FormerI].InnerEnd>Aligns[LatterI].InnerPos+Endurance)
		{
			int Diff=int(Aligns[FormerI].InnerEnd)-int(Aligns[LatterI].InnerPos);
			Diff/=2;
			Aligns[FormerI].InnerEnd-=Diff;
			Aligns[FormerI].End-=Diff;
			Aligns[FormerI].InnerLength-=Diff;
			Aligns[FormerI].Length-=Diff;
			Aligns[LatterI].InnerPos+=Diff;
			Aligns[LatterI].Pos+=Diff;
			Aligns[LatterI].InnerLength-=Diff;
			Aligns[LatterI].Length-=Diff;
		}
	}
}

void searchForClipSignatures(bam1_t *br, Contig & TheContig, Sam &SamFile, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, double *CoverageWindows, HandleBrMutex *Mut,Arguments & Args)
{
	if (!align_is_primary(br)) return;
	// int AS=bam_aux2i(bam_aux_get(br,"AS"));
	// int NM=bam_aux2i(bam_aux_get(br,"NM"));
	// int qlen=bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br));
	// int qlenc=brGetClippedQlen(br);//no clip qlen
	// if (bam_cigar_op(bam_get_cigar(br)[0])==BAM_CSOFT_CLIP) qlenc-=bam_cigar_oplen(bam_get_cigar(br)[0]);
	// if (bam_cigar_op(bam_get_cigar(br)[br->core.n_cigar-1])==BAM_CSOFT_CLIP) qlenc-=bam_cigar_oplen(bam_get_cigar(br)[br->core.n_cigar-1]);
	// printf("%d %d %d %d %d\n",AS,NM,qlenc,qlen,bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br)));
	// return;
	if (Tech==0) statCoverageCigar(br,CoverageWindows,TheContig,Mut,Args);
	vector<Alignment> Aligns;
	char* SA_tag_char = bam_get_string_tag(br, "SA");
	if(SA_tag_char == NULL)
	{
		// statCoverageCigar(br,CoverageWindows,TheContig,Args);
		return;
	}
	// printf("name:%s; SA:%s;\n",bam_get_qname(br),SA_tag_char);
	int SACharLen=strlen(SA_tag_char);
	char SATag[SACharLen+1];
	strcpy(SATag,SA_tag_char);
	for (int i=0;i<SACharLen;++i) if (SATag[i]==',' || SATag[i]==';') SATag[i]='\0';
	int k=0;
	int pos,strand,cid,mapQ,NM;
	const char * cigars;
	pos=br->core.pos;
	strand= read_is_forward(br)?1:0;
	Aligns.push_back(Alignment(br->core.tid,pos,strand,bam_get_cigar(br),br->core.n_cigar,br->core.qual,bam_aux2i(bam_aux_get(br,"NM"))));
	for (int i=0;i<SACharLen;i+=strlen(SATag+i)+1)
	{
		if (k%6==0)//rname
		{
			if (Args.CallByContig && strcmp(SATag+i,SamFile.Header->target_name[br->core.tid])!=0)//Pass other contig. Should be altered if want to do multi-chromosome sv.
			{
				++k;
				// i+=strlen(SATag+i)+1;
				for (;k%6!=0;++k) i+=strlen(SATag+i)+1;
				continue;
			}
			cid=Args.ContigNameID[SATag+i];
		}
		else if (k%6==1) pos=atoi(SATag+i);//pos
		else if (k%6==2) strand=SATag[i]=='+'?1:0;//strand
		else if (k%6==3) cigars=SATag+i;//CIGAR
		else if (k%6==4) mapQ=atoi(SATag+i);//mapQ
		else if (k%6==5)
		{
			NM=atoi(SATag+i);
			// printf(" %s",cigars);
			// if (abs(pos-br->core.pos)<=1000000)
			Aligns.push_back(Alignment(cid,pos,strand,cigars,mapQ,NM));
		}
		k+=1;
	}
	sort(Aligns.data(),Aligns.data()+Aligns.size());
	// dealClipConflicts(Aligns,Args);//Careful use. Good for ont but not for sims.
	// printf("name:%s; size:%lu;",bam_get_qname(br),Aligns.size());
	// printf("name:%s; %d,%d,%d;",bam_get_qname(br),Aligns[0].Strand,Aligns[0].Pos,Aligns[0].Length);
	// for (int i=0;i<Aligns.size();++i) printf(" %d,%d,%d",Aligns[i].Strand,Aligns[i].Pos,Aligns[i].Length);
	// for (int i=0;i<Aligns.size();++i) printf(" %d,%d,%d",Aligns[i].Strand,Aligns[i].ForwardPos,Aligns[i].ForwardEnd);
	// printf("\n");
	// for (int i=0;i<Aligns.size();++i)
	// {
	// 	statCoverage(Aligns[i].Pos, Aligns[i].End, CoverageWindows, TheContig, Args);//Still no other contigs.
	// 	// getDelFromCigar(Aligns[i].CIGAR.data(), Aligns[i].CIGAR.size(),Aligns[i].Pos, bam_get_qname(br), Tech, TypeSignatures[0], Args);
	// }
	searchDelFromAligns(br,TheContig,Aligns,Tech,TypeSignatures, Mut, Args);
	searchInsFromAligns(br,TheContig,Aligns,Tech,TypeSignatures, Mut, Args);
	searchDupFromAligns(br,TheContig,Aligns,Tech,TypeSignatures, Mut, Args);
	searchInvFromAligns(br,TheContig,Aligns,Tech,TypeSignatures, Mut, Args);
}

//This kind of signature should - some normal isize when calc svlen
void getDRPSignature(bam1_t * br, Stats& SampleStats, HandleBrMutex *Mut, Contig& TheContig, vector<vector<vector<Signature>>> &TypeSignatures)
{
	// const char * qname=bam_get_qname(br);
	// if (strcmp(qname, "HISEQ1:93:H2YHMBCXX:2:1112:15493:51859")==0)
	// {
	// 	string CIGAR=cigar2string(bam_get_cigar(br),br->core.n_cigar);
	// 	int a=1;
	// }
	if ((!read_is_unmapped(br)) && (!read_mate_is_unmapped(br)) && br->core.mtid==br->core.tid && br->core.isize>=0)
	{
		int ForwardBase=-1, ReverseBase, ForwardEnd, ReverseEnd, DupStart, DupEnd, DupSize;// ForwardBase|-->ForwardEnd  ReverseEnd<--|ReverseBase
		// if (read_is_read1(br) && read_is_forward(br) && (!read_mate_is_forward(br)) || !(read_is_read1(br)) && (!read_is_forward(br)) && read_mate_is_forward(br))
		if (read_is_forward(br) && (!read_mate_is_forward(br)))
		{
			ForwardBase=br->core.pos;
			ForwardEnd=br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br));
			ReverseBase=ForwardBase+br->core.isize;
			ReverseEnd=br->core.mpos;
			DupStart=ForwardBase;
			DupEnd=ReverseBase;
			DupSize=SampleStats.Mean-(DupEnd-DupStart);
			if (DupSize<0) DupSize=0;
			// uint32_t * cigars=bam_get_cigar(br);
			// int isize=br->core.mpos+bam_cigar2qlen(br->core.n_cigar,cigars)-int(br->core.pos);
			// if (br->core.isize>SampleStats.Mean+3*SampleStats.SD) TypeSignatures[0].push_back(Signature(1,1,0,br->core.pos,br->core.pos+br->core.isize,bam_get_qname(br),br->core.qual,Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize-SampleStats.Mean));
			// else if (br->core.isize<SampleStats.Mean-3*SampleStats.SD && br->core.isize>-1000) TypeSignatures[2].push_back(Signature(1,1,2,br->core.mpos,bam_endpos(br),bam_get_qname(br),br->core.qual,Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize>0?SampleStats.Mean-br->core.isize:abs(br->core.isize)+brGetClippedQlen(br)));
			// if (isize>SampleStats.Mean+3*SampleStats.SD) TypeSignatures[0].push_back(Signature(1,1,0,br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br)),br->core.mpos,bam_get_qname(br),br->core.qual,Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize-SampleStats.Mean));
			// else if (isize<SampleStats.Mean-3*SampleStats.SD && br->core.isize>-1000) TypeSignatures[2].push_back(Signature(1,1,2,br->core.mpos,bam_endpos(br),bam_get_qname(br),br->core.qual,Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize>0?SampleStats.Mean-br->core.isize:abs(br->core.isize)+brGetClippedQlen(br)));
		}
		else if ((!read_is_forward(br)) && read_mate_is_forward(br))
		{
			ForwardBase=br->core.mpos;
			ForwardEnd=br->core.pos+br->core.isize;
			ReverseBase=br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br));
			ReverseEnd=br->core.pos;
			DupStart=ReverseEnd;
			DupEnd=ForwardEnd;
			DupSize=DupEnd-DupStart;
		}
		if (ForwardBase!=-1)
		{
			int BaseGap=ReverseBase-ForwardBase;
			// if (BaseGap>100000)
			// {
			// 	string CIGAR=cigar2string(bam_get_cigar(br),br->core.n_cigar);
			// 	char * qname=bam_get_qname(br);
			// 	int a=1;
			// }
			int DelLength=BaseGap-SampleStats.Mean;
			if (BaseGap>SampleStats.Mean+5*SampleStats.SD)
			{
				pthread_mutex_lock(&Mut->m_Sig[0]);
				TypeSignatures[TheContig.ID][0].push_back(Signature(1,1,0,ForwardBase+int(SampleStats.Mean/2),ForwardBase+int(SampleStats.Mean/2)+DelLength,bam_get_qname(br),br->core.qual,Segment(ForwardBase,ForwardEnd),Segment(ReverseEnd,ReverseBase),DelLength));
				pthread_mutex_unlock(&Mut->m_Sig[0]);
			}
			else if (BaseGap<SampleStats.Mean-3*SampleStats.SD)
			{
				pthread_mutex_lock(&Mut->m_Sig[2]);
				TypeSignatures[TheContig.ID][2].push_back(Signature(1,1,2,DupStart,DupEnd,bam_get_qname(br),br->core.qual,Segment(ForwardBase,ForwardEnd),Segment(ReverseEnd,ReverseBase),DupSize));
				pthread_mutex_unlock(&Mut->m_Sig[2]);
			}
		}
	}
}
const char * BamBases="NACNGNNNTNNNNNNN";
void getInsFromCigar(bam1_t *br, int Tech, vector<Signature>& Signatures, HandleBrMutex *Mut, Arguments & Args)
{
	if (br->core.qual<Args.MinMappingQuality) return;
	int TLength= bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br));
	if (TLength<Args.MinTemplateLength) return;
	uint32_t * cigars=bam_get_cigar(br);
	int CurrentStart=-1, CurrentLength=0;
	int Begin=br->core.pos;
	int QueryBegin=0;
	int CurrentQueryStart=QueryBegin;
	string Allele;
	//int MergeDis=500;
	int MinMaxMergeDis=Args.InsMinMaxMergeDis;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	float MaxMergeDisPortion=Args.DelMaxMergePortion;
	for (int i=0;i<br->core.n_cigar;++i)
	{
		if (bam_cigar_op(cigars[i])==BAM_CINS && bam_cigar_oplen(cigars[i])>=Args.MinSVLen)
		{
			// int rlen=bam_cigar2rlen(1,cigars+i);
			int qlen=bam_cigar_oplen(cigars[i]);
			// printf("%d %d %s\n",Begin,rlen,bam_get_qname(br));
			if (CurrentStart==-1)
			{
				CurrentStart=Begin;
				CurrentLength=qlen;
				CurrentQueryStart=QueryBegin;
				Allele="";
			}
			else
			{
				if (Begin-CurrentStart-CurrentLength>=(CurrentLength*MaxMergeDisPortion>MinMaxMergeDis?CurrentLength*MaxMergeDisPortion:MinMaxMergeDis))
				{
					if(CurrentLength>=Args.MinSVLen)
					{
						Signature TempSignature(0,Tech,1,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br),br->core.qual,Allele.c_str());
						pthread_mutex_lock(&Mut->m_Sig[1]);
						Signatures.push_back(TempSignature);
						pthread_mutex_unlock(&Mut->m_Sig[1]);
					}
					CurrentStart=Begin;
					CurrentLength=qlen;
					CurrentQueryStart=QueryBegin;
					Allele="";
				}
				else
				{
					CurrentLength+=qlen;
				}
			}
			for (int j=0;j<qlen;++j)
			{
				Allele+=BamBases[bam_seqi(bam_get_seq(br),QueryBegin+j)];
			}
		}
		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Begin+=bam_cigar_oplen(cigars[i]);
		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==1||bam_cigar_op(cigars[i])==4||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) QueryBegin+=bam_cigar_oplen(cigars[i]);
		//Begin+=bam_cigar2rlen(1,cigars+i);
	}
	if (CurrentStart!=-1)
	{
		if(CurrentLength>=Args.MinSVLen)
		{
			Signature TempSignature(0,Tech,1,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br),br->core.qual,Allele.c_str());
			pthread_mutex_lock(&Mut->m_Sig[1]);
			Signatures.push_back(TempSignature);
			pthread_mutex_unlock(&Mut->m_Sig[1]);
		}
	}
}

inline void getDelFromCigar(bam1_t * br, int Tech, vector<Signature>& Signatures, HandleBrMutex *Mut, Arguments & Args)
{
	uint32_t * cigars=bam_get_cigar(br);
	unsigned n_cigar=br->core.n_cigar, pos=br->core.pos;
	const char * qname=bam_get_qname(br);
	int qual=br->core.qual;
// 	return getDelFromCigar(bam_get_cigar(br), br->core.n_cigar, br->core.pos, bam_get_qname(br), br->core.qual, Tech, Signatures, Mut, Args);
	int TLength= bam_cigar2qlen(n_cigar,cigars);
	if (TLength<Args.MinTemplateLength) return;
	int CurrentStart=-1, CurrentLength=0;
	int Begin=pos;
	int LastBegin=Begin;
	//int MergeDis=500;
	int MinMaxMergeDis=Args.DelMinMaxMergeDis;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	float MaxMergeDisPortion=Args.DelMaxMergePortion;
	double MergeScore=0;
	vector<Segment> ShortSigs;
	for (int i=0;i<n_cigar;++i)
	{
		if (bam_cigar_op(cigars[i])==BAM_CDEL && bam_cigar_oplen(cigars[i])>=Args.MinSVLen)
		{
			// int rlen=bam_cigar2rlen(1,cigars+i);
			int rlen=bam_cigar_oplen(cigars[i]);
			// printf("%d %d %s\n",Begin,rlen,bam_get_qname(br));
			if (CurrentStart==-1)
			{
				CurrentStart=Begin;
				CurrentLength=rlen;
				MergeScore=0;
			}
			else
			{
			// if (CurrentLength<500)
			// {
			// 	MinMaxMergeDis=50;
			// 	MaxMergeDisPortion=0;
			// }
			// else
			// {
			// 	MinMaxMergeDis=Args.DelMinMaxMergeDis;
			// 	MaxMergeDisPortion=Args.DelMaxMergePortion;
			// }
			if (Begin-CurrentStart-CurrentLength>=(CurrentLength*MaxMergeDisPortion>MinMaxMergeDis?CurrentLength*MaxMergeDisPortion:MinMaxMergeDis))
			{
				MergeScore=100-MergeScore;
				if(CurrentLength>=Args.MinSVLen)
				{
					Signature Temp(0,Tech,0,CurrentStart,CurrentStart+CurrentLength,qname,MergeScore);
					pthread_mutex_lock(&Mut->m_Sig[0]);
					Signatures.push_back(Temp);
					pthread_mutex_unlock(&Mut->m_Sig[0]);
				}
				// printf("%d %d %s\n",CurrentStart,CurrentLength,bam_get_qname(br));
				CurrentStart=Begin;
				CurrentLength=rlen;
				MergeScore=0;
			}
			else
			{
				CurrentLength+=rlen;
				MergeScore+=1;
			}
			}
		}
		else if (bam_cigar_op(cigars[i])==BAM_CDEL)
		{
			ShortSigs.push_back(Segment(Begin, bam_cigar_oplen(cigars[i])));
		}
		LastBegin=Begin;
		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Begin+=bam_cigar_oplen(cigars[i]);
		//Begin+=bam_cigar2rlen(1,cigars+i);
	}
	if (CurrentStart!=-1)
	{
		MergeScore=100-MergeScore*1;
		int AS=bam_aux2i(bam_aux_get(br,"AS"));
		int ClippedQLen=brGetClippedQlen(br);
		double Quality=((double)AS)/((double)ClippedQLen);
		if(CurrentLength>=Args.MinSVLen)
		{
			Signature TempSignature(0,Tech,0,CurrentStart,CurrentStart+CurrentLength,qname,Quality);
			pthread_mutex_lock(&Mut->m_Sig[0]);
			Signatures.push_back(TempSignature);
			pthread_mutex_unlock(&Mut->m_Sig[0]);
		}
	}
	//Merge short sigs
	// double MergeRatio=0.5;
	// int SkipGap=1000;
	// for (int i=0;i<ShortSigs.size();++i)
	// {
	// 	int Farthest=i;
	// 	int SigLength=ShortSigs[i].End-ShortSigs[i].Begin;
	// 	int TotalLength=0;
	// 	int Skipped=i;
	// 	for (int j=i+1;j<ShortSigs.size();++j)
	// 	{
	// 		SigLength+=ShortSigs[j].End-ShortSigs[j].Begin;
	// 		TotalLength=ShortSigs[j].End-ShortSigs[i].Begin;
	// 		if (double(SigLength)/double(TotalLength)>=MergeRatio)
	// 		{
	// 			Farthest=j;
	// 		}
	// 		if (ShortSigs[j].Begin>ShortSigs[j-1].End+SkipGap)
	// 		{
	// 			Skipped=j-1;
	// 			break;
	// 		}
	// 	}
	// 	if (Farthest!=i && ShortSigs[Farthest].End-ShortSigs[i].Begin>=Args.MinSVLen)
	// 	{
	// 		Signature TempSignature(0,Tech,0,ShortSigs[i].Begin,ShortSigs[Farthest].End,qname,100-Farthest+i);
	// 		pthread_mutex_lock(&Mut->m_Sig[0]);
	// 		Signatures.push_back(TempSignature);
	// 		pthread_mutex_unlock(&Mut->m_Sig[0]);
	// 		i=Farthest;
	// 	}
	// 	// if (Skipped!=i) i=Skipped;
	// }
}

// void getDelFromCigar(bam1_t *br, int Tech, vector<Signature>& Signatures, HandleBrMutex *Mut, Arguments & Args)
// {
// 	if (br->core.qual<Args.MinMappingQuality) return;
// 	return getDelFromCigar(bam_get_cigar(br), br->core.n_cigar, br->core.pos, bam_get_qname(br), br->core.qual, Tech, Signatures, Mut, Args);
// 	int TLength= bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br));
// 	if (TLength<Args.MinTemplateLength) return;
// 	uint32_t * cigars=bam_get_cigar(br);
// 	int CurrentStart=-1, CurrentLength=0;
// 	int Begin=br->core.pos;
// 	//int MergeDis=500;
// 	int MinMaxMergeDis=Args.DelMinMaxMergeDis;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
// 	float MaxMergeDisPortion=Args.DelMaxMergePortion;
// 	for (int i=0;i<br->core.n_cigar;++i)
// 	{
// 		if (bam_cigar_op(cigars[i])==BAM_CDEL && bam_cigar_oplen(cigars[i])>=Args.MinSVLen)
// 		{
// 			// int rlen=bam_cigar2rlen(1,cigars+i);
// 			int rlen=bam_cigar_oplen(cigars[i]);
// 			// printf("%d %d %s\n",Begin,rlen,bam_get_qname(br));
// 			if (CurrentStart==-1)
// 			{
// 				CurrentStart=Begin;
// 				CurrentLength=rlen;
// 			}
// 			else
// 			{
// 			if (Begin-CurrentStart-CurrentLength>=MIN(Args.DelMaxMaxMergeDis,MAX(CurrentLength*MaxMergeDisPortion,MinMaxMergeDis)))
// 			{
// 				if(CurrentLength>=Args.MinSVLen) Signatures.push_back(Signature(0,Tech,0,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br),br->core.qual));
// 				// printf("%d %d %s\n",CurrentStart,CurrentLength,bam_get_qname(br));
// 				CurrentStart=Begin;
// 				CurrentLength=rlen;
// 			}
// 			else
// 			{
// 				CurrentLength+=rlen;
// 			}
// 			}
// 		}
// 		//Those 4 kinds of op add ref
// 		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Begin+=bam_cigar_oplen(cigars[i]);
// 		//Begin+=bam_cigar2rlen(1,cigars+i);
// 	}
// 	if (CurrentStart!=-1)
// 	{
// 		if(CurrentLength>=Args.MinSVLen) Signatures.push_back(Signature(0,Tech,0,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br),br->core.qual));
// 				// printf("%d %d %s\n",CurrentStart,CurrentLength,bam_get_qname(br));
// 	}
// }

struct HandleBrArgs
{
	bam1_t *br;
	Contig *pTheContig;
	Sam *pSamFile;
	int Tech;
	Stats *pSampleStats;
	vector<vector<vector<Signature>>> *pTypeSignatures;
	SegmentSet *pAllPrimarySegments;
	double* CoverageWindows;
	HandleBrMutex *mut;
	Arguments * pArgs;
};

void handlebr(bam1_t *br, Contig * pTheContig, Sam *pSamFile, int Tech, Stats *pSampleStats, vector<vector<vector<Signature>>> *pTypeSignatures, SegmentSet * pAllPrimarySegments, double* CoverageWindows, HandleBrMutex *mut,Arguments * pArgs)
{
	Contig & TheContig=*pTheContig;
	Sam & SamFile=*pSamFile;
	Stats & SampleStats=*pSampleStats;
	SegmentSet & AllPrimarySegments=*pAllPrimarySegments;
	Arguments & Args=*pArgs;
	if (align_is_primary(br))
	{
		hts_pos_t end=br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br));
		pthread_mutex_lock(&mut->m_AllPrimarySeg);
		AllPrimarySegments.add(br->core.pos,end);
		pthread_mutex_unlock(&mut->m_AllPrimarySeg);
	}
	// int OldDELN=TypeSignatures[0].size();
	// int OldINSN=TypeSignatures[1].size();
	getDelFromCigar(br,Tech,(*pTypeSignatures)[pTheContig->ID][0], mut, Args);
	getInsFromCigar(br,Tech,(*pTypeSignatures)[pTheContig->ID][1], mut, Args);
	// int NewDELN=TypeSignatures[0].size()-OldDELN;
	// int NewINSN=TypeSignatures[1].size()-OldINSN;
	// if (NewINSN>0 && NewDELN>0 && ((NewDELN+NewINSN)>(MAX(10,br->core.l_qseq*0.001))))
	// {
	// 	for (int i=0;i<NewDELN;++i) TypeSignatures[0].pop_back();
	// 	for (int i=0;i<NewINSN;++i) TypeSignatures[1].pop_back();
	// }
	if (Tech==1)
	{
		if (read_is_paired(br))
		{
			getDRPSignature(br, SampleStats, mut, TheContig, *pTypeSignatures);
		}
	}
	searchForClipSignatures(br, TheContig, SamFile, Tech, *pTypeSignatures, CoverageWindows, mut, Args);
}

void * handlebrWrapper(void * args)
{
	HandleBrArgs * A=(HandleBrArgs*)args;
	handlebr(A->br, A->pTheContig, A->pSamFile, A->Tech, A->pSampleStats, A->pTypeSignatures, A->pAllPrimarySegments, A->CoverageWindows, A->mut, A->pArgs);
	bam_destroy1(A->br);
	delete A;
	return NULL;
}

/*
void handleBrBlock(BrBlock * TheBlock, Contig * Contigs, int * CordinTrans, int& ReadCount, int & UnmappedCount)
{
	int BlockReadCount=0, BlockUnmappedCount=0;
	int BrReadCount=0,BrUnmappedCount=0;
	for (int i=0;i<TheBlock->Count;++i)
	{
		BrReadCount=0;
		BrUnmappedCount=0;
		handlebr(TheBlock->Block+i,Contigs,CordinTrans,BrReadCount,BrUnmappedCount);
		BlockReadCount+=BrReadCount;
		BlockUnmappedCount+=BrUnmappedCount;
	}
	WriteLock.lock();
	ReadCount+=BlockReadCount;
	UnmappedCount+=BlockUnmappedCount;
	WriteLock.unlock();
}

void getBlockAndProcess(BrBlock ** UpMostBlock, Contig * Contigs, int * CordinTrans, int& ReadCount, int & UnmappedCount)
{
	while (true)
	{
		if ((*UpMostBlock)->Over) break;
		TopLock.lock();
		if (*UpMostBlock==NULL || !((*UpMostBlock)->Mature))
		{
			TopLock.unlock();
			sleep(1);
		}
		else
		{
			BrBlock * TheBlock=*UpMostBlock;
			if (TheBlock->Over)
			{
				TopLock.unlock();
				break;
			}
			else *UpMostBlock=(*UpMostBlock)->Next;
			TopLock.unlock();
			handleBrBlock(TheBlock,Contigs,CordinTrans,ReadCount,UnmappedCount);
			delete TheBlock;
		}
	}
}

void readBamToBrBlock(htsFile * SamFile,bam_hdr_t *Header, BrBlock** Top)
{
	BrBlock * TheBlock=new BrBlock();
	TheBlock->Next=new BrBlock();
	while (sam_read1(SamFile, Header, TheBlock->Block+TheBlock->Count) >=0)
	{
		++TheBlock->Count;
		if (TheBlock->Count>=TheBlock->Size)
		{
			TheBlock->Mature=true;
			if (*Top==NULL)//The only one who could fill NULL top
			{
				*Top=TheBlock;
			}
			TheBlock=TheBlock->Next;
			TheBlock->Next=new BrBlock();
		}
	}
	TheBlock->Mature=true;
	TheBlock->Next->Over=true;
}
*/

void takePipeAndHandleBr(Contig &TheContig, Sam & SamFile, int Tech, Stats & SampleStats, vector<vector<vector<Signature>>> &TypeSignatures, SegmentSet & AllPrimarySegments, double * CoverageWindows, HandleBrMutex *Mut, Arguments & Args, FILE* Pipe=stdin)
{
    bam1_t *br=bam_init1();
	size_t linebuffersize=1024*0124;
	char * linebuffer=(char*) malloc(linebuffersize);
	StdinLock.lock();
	int length=getline(&linebuffer,&linebuffersize,Pipe);
	StdinLock.unlock();
	//unsigned long long ReadCount=0;
	while (length>=0)
	{
		kstring_t ks={length,linebuffersize,linebuffer};
		sam_parse1(&ks,SamFile.Header,br);
		handlebr(br,&TheContig,&SamFile,Tech, &SampleStats, &TypeSignatures, &AllPrimarySegments, CoverageWindows, Mut, &Args);
		StdinLock.lock();
		length=getline(&linebuffer,&linebuffersize,Pipe);
		StdinLock.unlock();
	}
	free(linebuffer);
	bam_destroy1(br);
}

void calcMeanSD(StatsManager & SM, const char * SampleFileName, float &Mean, float &SD)
{
	Mean=0,SD=0;
	const int Size=10000;
	for (int i=0;i<=Size;++i)
	{
		Mean+=SM.getInsertLen(SampleFileName,float(i)/float(Size));
	}
	Mean/=float(Size+1);
	for (int i=0;i<=Size;++i)
	{
		float V=SM.getInsertLen(SampleFileName,float(i)/float(Size));
		//fprintf(stderr,"%f\n",V);
		V-=Mean;
		V*=V;
		SD+=V;
	}
	SD/=float(Size+1);
	SD=pow(SD,0.5);
}

Stats getSampleStats(const char * ReferenceFileName, const char * SampleFileName,int Tech=1)
{
		Stats SampleStats={0,0,0,0,0};
		if (Tech!=1) return SampleStats;
		StatsManager SM(ReferenceFileName,"");
		SM.handleBamCramStats(SampleFileName);
		SampleStats.MedianIS=SM.getInsertLen(SampleFileName,0.5f);
		SampleStats.UpIS=SM.getInsertLen(SampleFileName,0.99f);
		SampleStats.BelowIS=SM.getInsertLen(SampleFileName,0.01f);
		calcMeanSD(SM, SampleFileName, SampleStats.Mean, SampleStats.SD);
		return SampleStats;
}

vector<Stats> getAllStats(const char * ReferenceFileName, const vector<const char *> & BamFileNames, vector<int> AllTechs)
{
	vector<Stats> AllStats;
	for (int k=0;k<BamFileNames.size();++k)
	{
		const char * SampleFileName=BamFileNames[k];
		AllStats.push_back(getSampleStats(ReferenceFileName,SampleFileName,AllTechs[k]));
	}
	return AllStats;
}

string getSampleName(bam_hdr_t* Header)
{
	char SampleName[1024];
	kstring_t ks=KS_INITIALIZE;
	ks_resize(&ks,1000);
	sam_hdr_find_line_pos(Header,"RG",0,&ks);
	int splitn;
	int * splitoffsets=ksplit(&ks,'\t',&splitn);
	for (int i=0;i<splitn;++i)
	{
		if (memcmp(ks.s+splitoffsets[i],"SM",2)==0)
		{
			int End=ks_len(&ks);
			if (i<splitn-1) End=splitoffsets[i+1];
			memcpy(SampleName,ks.s+splitoffsets[i]+3,End-splitoffsets[i]-3);
			SampleName[End-splitoffsets[i]-3]='\0';
		}
	}
	ks_free(&ks);
	return string(SampleName);
}
Sam::Sam():SamFile(NULL), Header(NULL), BamIndex(NULL){}
void Sam::close()
{
	if (SamFile!=NULL) hts_close(SamFile);
}
htsThreadPool p = {NULL, 0};
vector<Sam> initSam(Arguments & Args)
{
	vector<Sam> SamFiles;
	if (Args.ThreadN > 1) {
		if (!(p.pool = hts_tpool_init(Args.ThreadN))) {
			die("Error creating thread pool\n");
		}
		p.qsize=Args.ThreadN*2;
	}
	for (int i=0;i<Args.BamFileNames.size();++i)
	{
		SamFiles.push_back(Sam());
		const char * SampleFileName=Args.BamFileNames[i];
		SamFiles[i].SamFile = hts_open(SampleFileName, "rb");
		SamFiles[i].BamIndex= sam_index_load(SamFiles[i].SamFile,SampleFileName);
		if (SamFiles[i].BamIndex==NULL)
		{
			die("Please index the samfile(%s) first!",SampleFileName);
		}
		//set reference file
		const char * ReferenceFileName=Args.ReferenceFileName;
		if (NULL != ReferenceFileName)
		{
			char referenceFilenameIndex[128];
			strcpy(referenceFilenameIndex, ReferenceFileName);
			strcat(referenceFilenameIndex, ".fai");
			int ret = hts_set_fai_filename(SamFiles[i].SamFile, referenceFilenameIndex);
		}
		SamFiles[i].Header = sam_hdr_read(SamFiles[i].SamFile);
		if (Args.SampleName=="*")
		{
			string SampleName=getSampleName(SamFiles[i].Header);
			if (SampleName!="") Args.SampleName=SampleName;
		}
		if (p.pool)
		{
			hts_set_opt(SamFiles[i].SamFile,  HTS_OPT_THREAD_POOL, &p);
			// if (settings.out) hts_set_opt(settings.out, HTS_OPT_THREAD_POOL, &p);
		}
	}
	return SamFiles;
}

void closeSam(vector<Sam> &SamFiles)
{
	for (int i=0;i<SamFiles.size();++i) SamFiles[i].close();
	if (p.pool) hts_tpool_destroy(p.pool);
}

void collectSignatures(Contig &TheContig, vector<vector<vector<Signature>>> &TypeSignatures, SegmentSet & AllPrimarySegments, Arguments & Args, vector<Sam>& SamFiles, vector<Stats> AllStats, vector<int> AllTechs, double* CoverageWindows, const char * DataSource)
{
	const char * ReferenceFileName=Args.ReferenceFileName;
	const vector<const char *> & BamFileNames=Args.BamFileNames;
	for (int k=0;k<BamFileNames.size();++k)
	{
		const char * SampleFileName=BamFileNames[k];
		fprintf(stderr,"Reading %s for region %s...\n",SampleFileName,TheContig.Name.c_str());
		Stats &SampleStats=AllStats[k];

		string Region=TheContig.Name;

		int Tech=AllTechs[k];

		int ReadCount=0, UnmappedCount=0;
		FILE * DSFile=0;
		HandleBrMutex mut;
		if (DataSource!=0)
		{
			if (strcmp(DataSource,"-")==0) DSFile=stdin;
			else
			{
				char *cmd=(char*) malloc(102400);
				cmd[0]='\0';
				if (strcmp(DataSource,"samtools")==0)
				{
					snprintf(cmd,102400,"samtools view -@ %d -T %s %s %s",Args.ThreadN-ReadThreadN-1,ReferenceFileName,SampleFileName,TheContig.Name.c_str());
					//fprintf(stderr,cmd);
				}
				DSFile = popen(cmd, "r");
				free(cmd);
				if (!DSFile) throw runtime_error("popen() failed!");
			}
		}
		if (DataSource!=0)
		{
			takePipeAndHandleBr(TheContig, SamFiles[k], Tech, SampleStats, TypeSignatures,AllPrimarySegments,CoverageWindows,&mut,Args,DSFile);
		}
		else
		{
			if (Args.ThreadN==1)
			{
				bam1_t *br=bam_init1();
				hts_itr_t* RegionIter=sam_itr_querys(SamFiles[k].BamIndex,SamFiles[k].Header,Region.c_str());
				while(sam_itr_next(SamFiles[k].SamFile, RegionIter, br) >=0)//read record
				{
					handlebr(br,&TheContig, &SamFiles[k], Tech, &SampleStats, &TypeSignatures, &AllPrimarySegments, CoverageWindows,&mut, &Args);
				}
				bam_destroy1(br);
			}
			else
			{
				hts_tpool_process *HandlebrProcess=hts_tpool_process_init(p.pool,p.qsize,1);
				bam1_t *br=bam_init1();
				hts_itr_t* RegionIter=sam_itr_querys(SamFiles[k].BamIndex,SamFiles[k].Header,Region.c_str());
				while(sam_itr_next(SamFiles[k].SamFile, RegionIter, br) >=0)//read record
				{
					bam1_t *cbr=bam_dup1(br);
					HandleBrArgs *A=new HandleBrArgs{cbr,&TheContig, &SamFiles[k], Tech, &SampleStats, &TypeSignatures, &AllPrimarySegments, CoverageWindows,&mut, &Args};
					hts_tpool_dispatch(p.pool,HandlebrProcess,handlebrWrapper,(void *)A);
				}
				bam_destroy1(br);
				hts_tpool_process_destroy(HandlebrProcess);
			}
		}
		if (DataSource!=0 && strcmp(DataSource,"-")!=0) pclose(DSFile);
	}
}

Contig * getContigs(Arguments & Args, int& NSeq, int RDWindowSize)
{
	faidx_t * Ref=fai_load(Args.ReferenceFileName);
	NSeq=faidx_nseq(Ref);
	Contig * Contigs=(Contig*) malloc(sizeof(Contig)*(NSeq));
	for (int i=0;i<(NSeq);++i)
	{
		const char * ContigName=faidx_iseq(Ref,i);
		int SeqLen=faidx_seq_len(Ref,ContigName);
		new (Contigs+i) Contig(i,ContigName, SeqLen);
		Args.ContigNameID[ContigName]=i;
	}
	fai_destroy(Ref);
	return Contigs;
}