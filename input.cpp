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
#include "defines.h"

using namespace std;
using namespace cre;

const int ReadThreadN=1;//old mechanism that leaves this number of threads handling reads. not good. obsoleted. use htslib's thread mechanisms instead
int ThreadN=8;

float MedianInsertionSize=481.2;
float ISSD=36.4;
//bool DynamicDRP=false;
//int DynamicCount=1000;
int RDWindowSize=100;

mutex TopLock, WriteLock, StdinLock;

#define read_is_unmapped(b) (((b)->core.flag&BAM_FUNMAP) != 0)
#define read_mate_is_unmapped(b) (((b)->core.flag&BAM_FMUNMAP) != 0)
#define read_is_paired(b) (((b)->core.flag&BAM_FPAIRED) != 0)
#define read_is_read1(b) (((b)->core.flag&BAM_FREAD1) != 0)
#define read_is_read2(b) (((b)->core.flag&BAM_FREAD2) != 0)
#define align_is_primary(b) ((((b)->core.flag&BAM_FSECONDARY) == 0) && (((b)->core.flag&BAM_FSUPPLEMENTARY) == 0))
#define align_is_secondary(b) (((b)->core.flag&BAM_FSECONDARY) != 0)
#define align_is_supplementary(b) (((b)->core.flag&BAM_FSUPPLEMENTARY) != 0)
#define read_is_forward(b) (((b)->core.flag&BAM_FREVERSE) == 0)
#define read_mate_is_forward(b) (((b)->core.flag&BAM_FMREVERSE) == 0)

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
	if (bam_cigar_op(CIGARD[CIGARN])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[CIGARN])==BAM_CHARD_CLIP)
	{
		--N;
		Clipped+=bam_cigar_oplen(CIGARD[CIGARN]);
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
	int Pos;
	int Length;//Length in Reference
	int End;
	int Strand;//0: -, 1: +
	vector<uint32_t> CIGAR;
	int InnerPos;
	int InnerLength;
	int InnerEnd;
	int ForwardPos;//for inner sort
	int ForwardEnd;
	Alignment(int Pos, int Strand, uint32_t* CIGARD, uint32_t CIGARN)
	:Pos(Pos), Strand(Strand)
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
	Alignment(int Pos, int Strand, const char * CIGARS)
	{
		size_t CIGARN=size_t(strlen(CIGARS)/2);
		uint32_t * ca=(uint32_t*) malloc(CIGARN*sizeof(uint32_t));
		sam_parse_cigar(CIGARS,NULL,&ca,&CIGARN);
		construct(Pos,Strand,ca,CIGARN);
	}
	bool operator<(const Alignment &other) const
	{
		return this->ForwardPos<other.ForwardPos;
	}
};

void searchDelFromAligns(bam1_t *br,vector<Alignment> &Aligns, int Tech, vector<Signature> &Signatures, Arguments & Args)
{
	for (int i=1;i<Aligns.size();++i)
	{
		//if (Aligns[i-1].Strand!=Aligns[i].Strand) continue;
		if (Aligns[i-1].End<Aligns[i].Pos && Aligns[i].Pos-Aligns[i-1].End-(Aligns[i].InnerPos-Aligns[i-1].InnerEnd)>=Args.MinSVLen) Signatures.push_back(Signature(2,Tech,0,Aligns[i-1].End,Aligns[i].Pos,bam_get_qname(br)));
	}
}

void searchDupFromAligns(bam1_t *br,vector<Alignment> &Aligns, int Tech, vector<Signature> &Signatures, Arguments & Args)
{
	for (int i=1;i<Aligns.size();++i)
	{
		if (Aligns[i-1].Strand!=Aligns[i].Strand) continue;
		vector<Segment> v;
		if (Aligns[i].Pos<Aligns[i-1].End)
		{
			int Dup=0;
			if (Aligns[i-1].End-Aligns[i].Pos>Aligns[i].Length)
			{
				Dup=1;
			}
			else if (Aligns[i-1].End-Aligns[i].Pos>Aligns[i-1].Length)
			{
				Dup=1;
			}
			else
			{
				if (Aligns[i-1].End-Aligns[i].Pos>=Args.MinSVLen)
				{
					Dup=1;
				}
			}
			if (Dup)
			{
				v.push_back(Segment(Aligns[i-1].Pos,Aligns[i-1].End));
				v.push_back(Segment(Aligns[i].Pos,Aligns[i].End));
				Signatures.push_back(Signature(2,Tech,1,Aligns[i].Pos,Aligns[i-1].End,bam_get_qname(br),v));
			}
		}
		//if (Aligns[i-1].End<Aligns[i].Pos && Aligns[i].Pos-Aligns[i-1].End-(Aligns[i].InnerPos-Aligns[i-1].InnerEnd)>=50) Signatures.push_back(Signature(2,Tech,1,Aligns[i-1].End,Aligns[i].Pos,bam_get_qname(br)));
	}
}

void searchForClipSignatures(bam1_t *br, Contig & TheContig, htsFile* SamFile, bam_hdr_t * Header, hts_idx_t* BamIndex, int Tech, vector<Signature> *TypeSignatures, Arguments & Args)
{
	if (!align_is_primary(br)) return;
	vector<Alignment> Aligns;
	char* SA_tag_char = bam_get_string_tag(br, "SA");
	if(SA_tag_char == NULL)
		return;
	int SACharLen=strlen(SA_tag_char)+1;
	char SATag[SACharLen];
	strcpy(SATag,SA_tag_char);
	for (int i=0;i<SACharLen;++i) if (SATag[i]==',' || SATag[i]==';') SATag[i]='\0';
	int k=0;
	int pos,strand;
	const char * cigars;
	pos=br->core.pos;
	strand= read_is_forward(br)?1:0;
	Aligns.push_back(Alignment(pos,strand,bam_get_cigar(br),br->core.n_cigar));
	for (int i=0;i<SACharLen;i+=strlen(SATag+i)+1)
	{
		if (k%6==0)//rname
		{
			if (strcmp(SATag+i,Header->target_name[br->core.tid])!=0)
			{
				++k;
				i+=strlen(SATag+i)+1;
				for (;k%6!=0;++k) i+=strlen(SATag+i)+1;
			}
		}
		else if (k%6==1) pos=atoi(SATag+i);//pos
		else if (k%6==2) strand=SATag[i]=='+'?1:0;//strand
		else if (k%6==3) cigars=SATag+i;//CIGAR
		else if (k%6==5)
		{
			if (abs(pos-br->core.pos)<=1000000)
				Aligns.push_back(Alignment(pos,strand,cigars));
		}
		k+=1;
	}
	sort(Aligns.data(),Aligns.data()+Aligns.size());
	searchDelFromAligns(br,Aligns,Tech,TypeSignatures[0], Args);
	searchDupFromAligns(br,Aligns,Tech,TypeSignatures[1], Args);
}

//This kind of signature should - some normal isize when calc svlen
void getDRPSignature(bam1_t * br, Stats& SampleStats, vector<Signature> *TypeSignatures)
{
	if ((!read_is_unmapped(br)) && (!read_mate_is_unmapped(br)))
	{
		if (read_is_forward(br) && (!read_mate_is_forward(br)))
		{
			//uint32_t * cigars=bam_get_cigar(br);
			if (br->core.isize>SampleStats.Mean+3*SampleStats.SD) TypeSignatures[0].push_back(Signature(1,1,0,br->core.pos,br->core.pos+br->core.isize,bam_get_qname(br),Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize-SampleStats.Mean));
			else if (br->core.isize<SampleStats.Mean-3*SampleStats.SD && br->core.isize>-1000) TypeSignatures[1].push_back(Signature(1,1,1,br->core.mpos,bam_endpos(br),bam_get_qname(br),Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize>0?SampleStats.Mean-br->core.isize:abs(br->core.isize)+brGetClippedQlen(br)));
		}
	}
}

void getDelFromCigar(bam1_t *br, int Tech, vector<Signature>& Signatures, Arguments & Args)
{
	#ifdef CUTE_VER
	if (br->core.qual<20) return;
	int TLength= bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br));
	if (TLength<500) return;
	uint32_t * cigars=bam_get_cigar(br);
	int CurrentStart=-1, CurrentLength=0;
	int Begin=br->core.pos;
	int MinMaxMergeDis=0;//Args.DelMinMaxMergeDis;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	float MaxMergeDisPortion=Args.DelMaxMergePortion;
	for (int i=0;i<br->core.n_cigar;++i)
	{
		if (bam_cigar_op(cigars[i])==BAM_CDEL && bam_cigar_oplen(cigars[i])>=Args.MinSVLen)
		{
			int rlen=bam_cigar_oplen(cigars[i]);
			if (CurrentStart==-1)
			{
				CurrentStart=Begin;
				CurrentLength=rlen;
			}
			else
			{
			if (Begin-CurrentStart-CurrentLength>=MinMaxMergeDis)//(CurrentLength*MaxMergeDisPortion>MinMaxMergeDis?CurrentLength*MaxMergeDisPortion:MinMaxMergeDis))
			{
				if(CurrentLength>=Args.MinSVLen) Signatures.push_back(Signature(0,Tech,0,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br)));
				CurrentStart=Begin;
				CurrentLength=rlen;
			}
			else
			{
				CurrentLength+=rlen;
			}
			}
		}
		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Begin+=bam_cigar_oplen(cigars[i]);
	}
	if (CurrentStart!=-1)
	{
		if(CurrentLength>=Args.MinSVLen) Signatures.push_back(Signature(0,Tech,0,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br)));
	}
	return;
	#endif
	// printf("%s %d %d\n", bam_get_qname(br),br->core.pos, br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br)));
	if (br->core.qual<20) return;
	int TLength= bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br));
	if (TLength<500) return;
	uint32_t * cigars=bam_get_cigar(br);
	int CurrentStart=-1, CurrentLength=0;
	int Begin=br->core.pos;
	//int MergeDis=500;
	int MinMaxMergeDis=0;//Args.DelMinMaxMergeDis;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	float MaxMergeDisPortion=Args.DelMaxMergePortion;
	for (int i=0;i<br->core.n_cigar;++i)
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
			}
			else
			{
			if (Begin-CurrentStart-CurrentLength>=MinMaxMergeDis)//(CurrentLength*MaxMergeDisPortion>MinMaxMergeDis?CurrentLength*MaxMergeDisPortion:MinMaxMergeDis))
			{
				if(CurrentLength>=Args.MinSVLen) Signatures.push_back(Signature(0,Tech,0,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br)));
				// printf("%d %d %s\n",CurrentStart,CurrentLength,bam_get_qname(br));
				CurrentStart=Begin;
				CurrentLength=rlen;
			}
			else
			{
				CurrentLength+=rlen;
			}
			}
		}
		// if (bam_cigar_op(cigars[i])==BAM_CINS)
		// {
		// 	if (CurrentStart!=-1)
		// 	{
		// 		int rlen=bam_cigar2rlen(1,cigars+i);
		// 		CurrentLength-=rlen;
		// 	}
		// }
		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Begin+=bam_cigar_oplen(cigars[i]);
		//Begin+=bam_cigar2rlen(1,cigars+i);
	}
	if (CurrentStart!=-1)
	{
		if(CurrentLength>=Args.MinSVLen) Signatures.push_back(Signature(0,Tech,0,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br)));
				// printf("%d %d %s\n",CurrentStart,CurrentLength,bam_get_qname(br));
	}
	/*
	vector<int> Splits;
	int Length=bam_cigar2rlen(br->core.n_cigar,cigars);
	if (Length<=0) return;
	int * Variants=(int*) calloc(Length, sizeof(int));
	int k=0;
	for (int i=0;i<br->core.n_cigar;++i)
	{
		int rlen=bam_cigar2rlen(1,cigars+i);
		int qlen=bam_cigar2qlen(1,cigars+i);
		for (int j=k;j<k+rlen;++j)
		{
			if (rlen<qlen) Variants[i]=1;
			else if (rlen>qlen) Variants[i]=-1;
			else Variants[i]=0;
		}
		k+=rlen;
	}
	int Window=5,ThresholdN=3,ThresholdSum=2, Last=0;
	for (int i=0;i<Length;++i)//find splits
	{
		int Sum=0,N=0;
		for (int j=i;j<i+Window;++j)
		{
			Sum+=Variants[j];
			if (Variants[j]!=0) N+=1;
			if (i+Window<Length)
			{
				if (N==Length-i && abs(Sum)>=1)
				{
					if (Last==0)
					{
						Splits.push_back(i);
						Last=1;
					}
				}
				else
				{
					if (Last==1)
					{
						Splits.push_back(i);
						Last=0;
					}
				}
			}
			else
			{
				if (N>=ThreadN && abs(Sum)>=ThresholdSum)
				{
					if (Last==0)
					{
						Splits.push_back(i);
						Last=1;
					}
				}
				else
				{
					if (Last==1)
					{
						Splits.push_back(i);
						Last=0;
					}
				}
			}
		}
		Last=0,Sum=0;
		for (int i=0;i<Splits.size();++i)
		{
			for (int j=Last;j<Splits[i];++j)
			{
				Sum+=Variants[j];
			}
			if (Sum<=-50) Signatures.push_back(Signature(0,Tech,0,br->core.pos+Last,br->core.pos+Splits[i],bam_get_qname(br)));
			Last=Splits[i];
		}
	}*/
		//if ((cigars[i]&0xf==BAM_CDEL)&&cigars[i]>>4>=50) Signatures.push_back(Signature(0,Tech,0,Begin,Begin+cigars[i]>>4,bam_get_qname(br)));
		//if (align_is_primary(br)&&((cigars[i]&0xf==BAM_CSOFT_CLIP)||cigars[i]&0xf==BAM_CHARD_CLIP)&&(cigars[i]>>4>=10)) searchForClipSignatures(br,TheContig, SamFile, Header, BamIndex, Tech, Signatures);
		//Begin+=bam_cigar2rlen(1,cigars+i);
}

void handlebr(bam1_t *br, Contig & TheContig, htsFile* SamFile, bam_hdr_t * Header, hts_idx_t* BamIndex, int Tech, Stats &SampleStats, vector<Signature> *TypeSignatures, Arguments & Args)
{
	#ifdef CUTE_VER
	getDelFromCigar(br,Tech,TypeSignatures[0],Args);
	return;
	#endif
	getDelFromCigar(br,Tech,TypeSignatures[0],Args);
	if (Tech==1)
	{
		if (read_is_paired(br) && read_is_read1(br))
		{
			getDRPSignature(br, SampleStats, TypeSignatures);
		}
	}
	// searchForClipSignatures(br, TheContig, SamFile, Header, BamIndex, Tech, TypeSignatures, Args);
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

void takePipeAndHandleBr(Contig &TheContig, htsFile* SamFile, bam_hdr_t * Header, hts_idx_t* BamIndex, int Tech, Stats & SampleStats, vector<Signature>* TypeSignatures, Arguments & Args, FILE* Pipe=stdin)
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
		sam_parse1(&ks,Header,br);
		handlebr(br,TheContig,SamFile,Header,BamIndex,Tech, SampleStats, TypeSignatures, Args);
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

void collectSignatures(Contig &TheContig, vector<Signature> *ContigTypeSignatures, Arguments & Args, vector<Stats> AllStats, vector<int> AllTechs, const char * DataSource)
{
	const char * ReferenceFileName=Args.ReferenceFileName;
	const vector<const char *> & BamFileNames=Args.BamFileNames;
	for (int k=0;k<BamFileNames.size();++k)
	{
		const char * SampleFileName=BamFileNames[k];
		fprintf(stderr,"Reading %s for region %s...\n",SampleFileName,TheContig.Name.c_str());
		Stats &SampleStats=AllStats[k];

		htsFile* SamFile;//the file of BAM/CRAM
		bam_hdr_t* Header;//header for BAM/CRAM file
		hts_idx_t* BamIndex;
		string Region=TheContig.Name;

		//if ((!filesystem::exists(string(SampleFileName)+".")) &&(!filesystem::exists(string(SampleFileName)+".")) &&(!filesystem::exists(string(SampleFileName)+".")))
		//{
		//	die("Please index the samfile(%s) first!",SampleFileName);
		//}
		SamFile = hts_open(SampleFileName, "rb");
		BamIndex= sam_index_load(SamFile,SampleFileName);
		if (BamIndex==NULL)
		{
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

		int Tech=AllTechs[k];
		if (Args.SampleName=="*")
		{
			string SampleName=getSampleName(Header);
			if (SampleName!="") Args.SampleName=SampleName;
		}

		int ReadCount=0, UnmappedCount=0;
		FILE * DSFile=0;
		if (DataSource!=0)
		{
			if (strcmp(DataSource,"-")==0) DSFile=stdin;
			else
			{
				char *cmd=(char*) malloc(102400);
				cmd[0]='\0';
				if (strcmp(DataSource,"samtools")==0)
				{
					snprintf(cmd,102400,"samtools view -@ %d -T %s %s %s",ThreadN-ReadThreadN-1,ReferenceFileName,SampleFileName,TheContig.Name.c_str());
					//fprintf(stderr,cmd);
				}
				DSFile = popen(cmd, "r");
				if (!DSFile) throw runtime_error("popen() failed!");
			}
		}
		if (DataSource!=0)
		{
			takePipeAndHandleBr(TheContig, SamFile, Header, BamIndex, Tech, SampleStats, ContigTypeSignatures,Args,DSFile);
		}
		else
		{
    		htsThreadPool p = {NULL, 0};
			if (ThreadN > 1) {
				if (!(p.pool = hts_tpool_init(ThreadN))) {
					die("Error creating thread pool\n");
				}
				hts_set_opt(SamFile,  HTS_OPT_THREAD_POOL, &p);
				// if (settings.out) hts_set_opt(settings.out, HTS_OPT_THREAD_POOL, &p);
			}
			bam1_t *br=bam_init1();
			hts_itr_t* RegionIter=sam_itr_querys(BamIndex,Header,Region.c_str());
			while(sam_itr_next(SamFile, RegionIter, br) >=0)//read record
			{
				handlebr(br,TheContig, SamFile, Header, BamIndex, Tech, SampleStats, ContigTypeSignatures, Args);
			}
			bam_destroy1(br);
			if (p.pool) hts_tpool_destroy(p.pool);
		}
		if (DataSource!=0 && strcmp(DataSource,"-")!=0) pclose(DSFile);
		hts_close(SamFile);
	}
}

Contig * getContigs(const char *ReferenceFileName, int& NSeq, int RDWindowSize)
{
	faidx_t * Ref=fai_load(ReferenceFileName);
	NSeq=faidx_nseq(Ref);
	Contig * Contigs=(Contig*) malloc(sizeof(Contig)*(NSeq));
	for (int i=0;i<(NSeq);++i)
	{
		const char * ContigName=faidx_iseq(Ref,i);
		int SeqLen=faidx_seq_len(Ref,ContigName);
		//int Size=SeqLen/RDWindowSize;
		//if (SeqLen%RDWindowSize!=0) ++Size;
		new (Contigs+i) Contig(ContigName, SeqLen);
	}
	fai_destroy(Ref);
	return Contigs;
}