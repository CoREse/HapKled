#include <vector>
#include "variant.h"
#include "signature.h"
#include "optutils/OptHelper.h"
#include "contig.h"
#include "input.h"
#include "kled.h"
#include "clustering.h"
#include "report.h"
#include "htslib/htslib/faidx.h"
#include <algorithm>
#include "crelib/crelib.h"
#include <omp.h>
#include <iterator>
#include <list>
#include "ThreadPool.h"
#include <unordered_map>
using namespace std;
using namespace cre;

void sortAndDeDup(vector<Signature> &V)
{
	if (V.size()==0) return;
	sort(V.begin(),V.begin()+V.size());
	for (int i=1;i<V.size();++i)
	{
		if (V[i]==V[i-1]) V[i].Type=-1;
	}
}

#pragma omp declare reduction(RecordVectorConc: vector<VCFRecord>: omp_out.insert(omp_out.end(),make_move_iterator(omp_in.begin()),make_move_iterator(omp_in.end())))
#pragma omp declare reduction(RecordListConc: list<VCFRecord>: omp_out.splice(omp_out.end(),omp_in))

Arguments Args;
int main(int argc, const char* argv[])
{
	bool NoHeader=false;
	OptHelper OH=OptHelper("kled [Options] Bam1 [Bam2] [Bam3] ...");
    OH.addOpt('N', 0, 1, "TestNumber", "for test notation",'i',&(Args.TestN));
    OH.addOpt('R', "Ref", 1, "FileName", "Indicate Reference Fasta File(required)",'s',&(Args.ReferenceFileName));
    OH.addOpt('C', 0, 1, "ContigName", "Only call variants in Contig(s), can occur multiple times",'s',&(Args.CallingContigs),true);
    OH.addOpt('S', 0, 1, "SampleName", "Sample name, if not given, kled will try to get it from the first bam file",'S',&(Args.SampleName));
    OH.addOpt(0, "NOH", 0, "", "No header, for test",'b',&(NoHeader));
    OH.addOpt(0, "CCS", 0, "", "All bams are CCS data.",'b',&(Args.AllCCS));
    OH.addOpt(0, "CLR", 0, "", "All bams are CLR data.",'b',&(Args.AllCLR));
    OH.addOpt(0, "NOF", 0, "", "No filter, output all results.(default false)",'b',&(Args.NoFilter));
    OH.getOpts(argc,argv);

	Args.BamFileNames=OH.Args;

	if (Args.ReferenceFileName==0 || Args.BamFileNames.size()==0)
	{
		OH.showhelp();
		exit(1);
	}

	omp_set_num_threads(8);
	// ThreadPool ThePool(8);

	updateTime("Starting kled, reading reference...");
	int NSeq;
	Contig * Contigs=getContigs(Args.ReferenceFileName,NSeq);//,RDWindowSize);
    
	updateTime("Reading reference","Getting stats...");
	vector<int> AllTechs=getAllTechs(Args);

	vector<Stats> AllStats=getAllStats(Args.ReferenceFileName,Args.BamFileNames,AllTechs);

	for (int i=0;i<AllStats.size();++i) fprintf(stderr,"%f %f %f %f %f\n",AllStats[i].BelowIS,AllStats[i].MedianIS,AllStats[i].UpIS,AllStats[i].Mean,AllStats[i].SD);
	//exit(0);

	VCFHeader Header(Args.ReferenceFileName);
	addKledEntries(Header);
	for (int i=0;i<NSeq;++i)
	{
		if (Args.CallingContigs.size()!=0)
		{
			bool ToCall=false;
			for (int k=0;k<Args.CallingContigs.size();++k)
			{
				if (Contigs[i].Name==Args.CallingContigs[k])
				{
					ToCall=true;
					break;
				}
			}
			if (!ToCall) continue;
		}
		Header.addContig(Contigs[i]);
	}

	faidx_t * Ref=fai_load(Args.ReferenceFileName);
	vector<vector<Variant>> VariantsByContig;
	bool FirstBam=true;
	vector<Sam> SamFiles=initSam(Args);
	double TotalCoverage=0;//Accumulative
	unsigned ProcessedLength=0;
	updateTime("Getting stats", "Starting calling...");
	for (int i=0;i<NSeq;++i)
	{
		if (Args.CallingContigs.size()!=0)
		{
			bool ToCall=false;
			for (int k=0;k<Args.CallingContigs.size();++k)
			{
				if (Contigs[i].Name==Args.CallingContigs[k])
				{
					ToCall=true;
					break;
				}
			}
			if (!ToCall) continue;
		}
		updateTime("","Calling...");
		vector<Signature> ContigTypeSignatures[NumberOfSVTypes];//For supported SV type
		unsigned int CoverageWindowSize=Args.CoverageWindowSize;
		unsigned int NumberOfCoverageWindows=Contigs[i].Size/CoverageWindowSize+1;
		double *CoverageWindows=new double[NumberOfCoverageWindows];
		for (int k=0;k<Contigs[i].Size/CoverageWindowSize+1;++k) CoverageWindows[k]=0;
		collectSignatures(Contigs[i],ContigTypeSignatures,Args,SamFiles,AllStats,AllTechs,CoverageWindows,0);
		fprintf(stderr,"%ld\n",Contigs[i].Size-1);
		double *CoverageWindowsSums=NULL;//=(double*) malloc(sizeof(double)*(int)(NumberOfCoverageWindows+1));
		CoverageWindows[0]=0;
		// CoverageWindowsSums[0]=0;
		int CheckPointInterval=10000;
		double *CheckPoints=NULL;//=(double *)malloc(sizeof(double)*(int)(NumberOfCoverageWindows/CheckPointInterval+1));
		// CheckPoints[0]=0;
		// for (int i=1;i<NumberOfCoverageWindows+1;++i)
		// {
		// 	CoverageWindowsSums[i]=CoverageWindowsSums[i-1]+CoverageWindows[i];
		// 	if (i%CheckPointInterval==0)
		// 	{
		// 		CheckPoints[(int)i/CheckPointInterval]=CoverageWindowsSums[i];
		// 		CoverageWindowsSums[i]=0;
		// 	}
		// }
		double WholeCoverage=getAverageCoverage(0,Contigs[i].Size-1,CoverageWindows,Args, CoverageWindowsSums, CheckPoints, CheckPointInterval);
		TotalCoverage=TotalCoverage*((double)(ProcessedLength)/(double)(ProcessedLength+Contigs[i].Size));
		TotalCoverage+=WholeCoverage*((double)(Contigs[i].Size)/(double)(ProcessedLength+Contigs[i].Size));
		ProcessedLength+=Contigs[i].Size;
		// double WholeCoverage=CoverageWindowsSums[(int)(Contigs[i].Size/CoverageWindowSize+1)]/(Contigs[i].Size/CoverageWindowSize+1);
		// continue;
		if (!NoHeader and FirstBam)
		{
			Header.addSample(Args.SampleName.c_str());
			printf(Header.genHeader().c_str());
			FirstBam=false;
		}
		int totalsig=0,cigardel=0, cigarins=0, cigardup=0, drpdel=0, drpdup=0, clipdel=0, clipins=0, clipdup=0, clipinv=0;
		for (int m=0;m<NumberOfSVTypes;++m)
		{
			vector<Signature>& ContigSignatures=ContigTypeSignatures[m];
			totalsig+=ContigSignatures.size();
			for (int j=0;j<ContigSignatures.size();++j)
			{
				if (ContigSignatures[j].Type==0)
				{
					if (ContigSignatures[j].SupportedSV==0) ++cigardel;
					if (ContigSignatures[j].SupportedSV==1) ++cigarins;
					if (ContigSignatures[j].SupportedSV==2) ++cigardup;
				}
				else if (ContigSignatures[j].Type==1)
				{
					if (ContigSignatures[j].SupportedSV==0) ++drpdel;
					if (ContigSignatures[j].SupportedSV==2) ++drpdup;
				}
				else
				{
					if (ContigSignatures[j].SupportedSV==0) ++clipdel;
					if (ContigSignatures[j].SupportedSV==1) ++clipins;
					if (ContigSignatures[j].SupportedSV==2) ++clipdup;
					if (ContigSignatures[j].SupportedSV==3) ++clipinv;
				}
			}
		}
		fprintf(stderr,"%s: %llu\n, cigardel: %d, cigarins: %d, cigardup: %d, drpdel: %d, drpdup: %d, clipdel: %d, clipins: %d, clipdup: %d, clipinv: %d. Contig Size:%ld, Average Coverage: %lf, Total Average Coverage: %lf\n",Contigs[i].Name.c_str(),totalsig,cigardel, cigarins, cigardup, drpdel, drpdup, clipdel, clipins, clipdup, clipinv, Contigs[i].Size, WholeCoverage, TotalCoverage);

		updateTime("Getting signatures","Clustering...");
		vector<vector<Signature>> SignatureTypeClusters[NumberOfSVTypes];
		for (int k=0;k<NumberOfSVTypes;++k)
		{
			sortAndDeDup(ContigTypeSignatures[k]);
			clustering(ContigTypeSignatures[k],SignatureTypeClusters[k],AllStats[i],Args);
		}
		updateTime("Clustering","Generating results...");
		vector<vector<Signature>> SignatureClusters;
		for (int k=0;k<NumberOfSVTypes;++k) SignatureClusters.insert(SignatureClusters.end(),make_move_iterator(SignatureTypeClusters[k].begin()),make_move_iterator(SignatureTypeClusters[k].end()));
		vector<VCFRecord> Records;
		// int Times[8]={0,0,0,0,0,0,0,0};
		// unordered_map<thread::id,vector<VCFRecord>> ThreadRecords;
		// list<future<vector<VCFRecord>>> ThreadResults;
		// int BatchSize=100;
		// for (int j=0;j<SignatureClusters.size();j+=BatchSize)
		// {
		// 	// ThreadResults.push_back(ThePool.enqueue([](vector<vector<Signature>> &SignatureClusters,int Start, int End,const Contig & TheContig, faidx_t * Ref, double* CoverageWindows, double WholeCoverage, Arguments& Args, double * CoverageWindowsSums, double * CheckPoints, int CheckPointInterval) {
		// 	ThreadResults.push_back(ThePool.enqueue([&]() {
		// 		// if (ThreadRecords.count(this_thread::get_id())==0) {ThreadRecords[this_thread::get_id()]=vector<VCFRecord>();
		// 		// fprintf(stderr, "%lu\n", ThreadRecords.hash_function()(this_thread::get_id()));}
		// 		vector<VCFRecord> Results;
		// 		for (int ti=j;ti<MIN(j+BatchSize,SignatureClusters.size());++ti)
		// 		{
		// 			Results.push_back(VCFRecord(Contigs[i],Ref,SignatureClusters[ti],CoverageWindows, WholeCoverage, Args, CoverageWindowsSums, CheckPoints, CheckPointInterval));
		// 		}
		// 		return Results;
		// 	}
		// 	// SignatureClusters, j, MIN(j+BatchSize,SignatureClusters.size()), Contigs[i], Ref, CoverageWindows, WholeCoverage, Args, CoverageWindowsSums, CheckPoints, CheckPointInterval
		// 	));
		// }
		// for (list<future<vector<VCFRecord>>>::iterator iter=ThreadResults.begin();iter!=ThreadResults.end();++iter)
		// {
		// 	vector<VCFRecord> Results=iter->get();
		// 	fprintf(stderr,"%d\n",Results.size());
		// 	Records.insert(Records.end(),make_move_iterator(Results.begin()),make_move_iterator(Results.end()));
		// }
		// for (unordered_map<thread::id,vector<VCFRecord>>::iterator iter=ThreadRecords.begin();iter!=ThreadRecords.end();++iter)
		// {
		// 	fprintf(stderr,"%d\n",iter->second.size());
		// 	Records.insert(Records.end(),make_move_iterator(iter->second.begin()),make_move_iterator(iter->second.end()));
		// }
		#pragma omp parallel for reduction(RecordVectorConc:Records)
		for (int j=0;j<SignatureClusters.size();++j)
		{
			// ++Times[omp_get_thread_num()];
			Records.push_back(VCFRecord(Contigs[i],Ref,SignatureClusters[j],CoverageWindows, TotalCoverage, Args, CoverageWindowsSums, CheckPoints, CheckPointInterval));
		}
		// fprintf(stderr,"Times: %d %d %d %d %d %d %d %d\n",Times[0],Times[1],Times[2],Times[3],Times[4],Times[5],Times[6],Times[7]);
		// for (vector<VCFRecord>::iterator iter=Records.begin();iter!=Records.end();++iter)
		// {
		// 	if (iter->Keep)
		// 	{
		// 		KeptRecords.insert(KeptRecords.end(),make_move_iterator(iter),make_move_iterator(iter+1));
		// 	}
		// }
		updateTime("Results generation","Sorting results...");
		sort(Records.data(),Records.data()+Records.size());
		// Records.sort();
		updateTime("Results sorting","Outputing results...");
		for (auto r: Records)
		{
			if (!r.Keep) continue;
			r.resolveRef(Contigs[i],Ref);
			printf("\n%s",string(r).c_str());
		}
		//vector<Variant> ContigVariants;
		//VariantsByContig.push_back(ContigVariants);
		//callVariants(Contigs[i],VariantsByContig[VariantsByContig.size()-1],ContigSignatures,Args);
		delete CoverageWindows;
		// free(CoverageWindowsSums);
		// free(CheckPoints);
	}

	//report(VariantsByContig);
	fai_destroy(Ref);
	free(Contigs);
	closeSam(SamFiles);
    return 0;
}