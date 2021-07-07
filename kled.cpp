#include <vector>
#include "variant.h"
#include "signature.h"
#include "optutils/OptHelper.h"
#include "contig.h"
#include "input.h"
#include "kled.h"
using namespace std;

Arguments Args;
int main(int argc, const char* argv[])
{
	OptHelper OH=OptHelper("kled [Options] Bam1 [Bam2] [Bam3] ...");
    OH.addOpt('R', "Ref", 1, "FileName", "Indicate Reference Fasta File",'s',&(Args.ReferenceFileName));
    OH.getOpts(argc,argv);

	Args.BamFileNames=OH.Args;

	if (Args.ReferenceFileName==0 || Args.BamFileNames.size()==0)
	{
		OH.showhelp();
		exit(1);
	}

	int NSeq;
	Contig * Contigs=getContigs(Args.ReferenceFileName,NSeq);//,RDWindowSize);
    
	vector<Stats> AllStats=getAllStats(Args.ReferenceFileName,Args.BamFileNames);

	for (int i=0;i<AllStats.size();++i) fprintf(stderr,"%f %f %f %f %f\n",AllStats[i].BelowIS,AllStats[i].MedianIS,AllStats[i].UpIS,AllStats[i].Mean,AllStats[i].SD);
	//exit(0);

	vector<vector<Variant>> VariantsByContig;
	for (int i=0;i<NSeq;++i)
	{
		vector<Signature> ContigSignatures;
		collectSignatures(Contigs[i],ContigSignatures,Args.ReferenceFileName,Args.BamFileNames,AllStats);
		int cigardel=0, cigardup=0, drpdel=0, drpdup=0, clipdel=0, clipdup=0;
		for (int i=0;i<ContigSignatures.size();++i)
		{
			if (ContigSignatures[i].Type==0)
			{
				if (ContigSignatures[i].SupportedSV==0) ++cigardel;
				if (ContigSignatures[i].SupportedSV==1) ++cigardup;
			}
			else if (ContigSignatures[i].Type==1)
			{
				if (ContigSignatures[i].SupportedSV==0) ++drpdel;
				if (ContigSignatures[i].SupportedSV==1) ++drpdup;
			}
			else
			{
				if (ContigSignatures[i].SupportedSV==0) ++clipdel;
				if (ContigSignatures[i].SupportedSV==1) ++clipdup;
			}
		}
		fprintf(stderr,"%s: %llu\n, cigardel: %d, cigardup: %d, drpdel: %d, drpdup: %d, clipdel: %d, clipdup: %d\n",Contigs[i].Name.c_str(),ContigSignatures.size(),cigardel, cigardup, drpdel, drpdup, clipdel, clipdup);
		vector<Variant> ContigVariants;
		VariantsByContig.push_back(ContigVariants);
		//callVariants(Contigs[i],VariantsByContig[VariantsByContig.size()-1],ContigSignatures,Args);
	}

	//report(VariantsByContig);
	
	free(Contigs);
    return 0;
}