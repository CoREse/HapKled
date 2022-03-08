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
using namespace std;

void sortAndDeDup(vector<Signature> &V)
{
	if (V.size()==0) return;
	sort(V.begin(),V.begin()+V.size());
	for (int i=1;i<V.size();++i)
	{
		if (V[i]==V[i-1]) V[i].Type=-1;
	}
}

Arguments Args;
int main(int argc, const char* argv[])
{
	bool NoHeader=false;
	OptHelper OH=OptHelper("kled [Options] Bam1 [Bam2] [Bam3] ...");
    OH.addOpt('R', "Ref", 1, "FileName", "Indicate Reference Fasta File(required)",'s',&(Args.ReferenceFileName));
    OH.addOpt('C', 0, 1, "ContigName", "Only call variants in Contig(s), can occur multiple times",'s',&(Args.CallingContigs),true);
    OH.addOpt('S', 0, 1, "SampleName", "Sample name, if not given, kled will try to get it from the first bam file",'S',&(Args.SampleName));
    OH.addOpt(0, "NOH", 0, "SampleName", "No header, for test",'b',&(NoHeader));
    OH.getOpts(argc,argv);

	Args.BamFileNames=OH.Args;

	if (Args.ReferenceFileName==0 || Args.BamFileNames.size()==0)
	{
		OH.showhelp();
		exit(1);
	}

	int NSeq;
	Contig * Contigs=getContigs(Args.ReferenceFileName,NSeq);//,RDWindowSize);
    
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
		vector<Signature> ContigTypeSignatures[2];//For supported SV type
		collectSignatures(Contigs[i],ContigTypeSignatures,Args,AllStats,AllTechs,0);
		// continue;
		if (!NoHeader and FirstBam)
		{
			Header.addSample(Args.SampleName.c_str());
			printf(Header.genHeader().c_str());
			FirstBam=false;
		}
		int cigardel=0, cigardup=0, drpdel=0, drpdup=0, clipdel=0, clipdup=0;
		for (int m=0;m<2;++m)
		{
			vector<Signature>& ContigSignatures=ContigTypeSignatures[m];
			for (int j=0;j<ContigSignatures.size();++j)
			{
				if (ContigSignatures[j].Type==0)
				{
					if (ContigSignatures[j].SupportedSV==0) ++cigardel;
					if (ContigSignatures[j].SupportedSV==1) ++cigardup;
				}
				else if (ContigSignatures[j].Type==1)
				{
					if (ContigSignatures[j].SupportedSV==0) ++drpdel;
					if (ContigSignatures[j].SupportedSV==1) ++drpdup;
				}
				else
				{
					if (ContigSignatures[j].SupportedSV==0) ++clipdel;
					if (ContigSignatures[j].SupportedSV==1) ++clipdup;
				}
			}
		}
		fprintf(stderr,"%s: %llu\n, cigardel: %d, cigardup: %d, drpdel: %d, drpdup: %d, clipdel: %d, clipdup: %d\n",Contigs[i].Name.c_str(),ContigTypeSignatures[0].size()+ContigTypeSignatures[1].size(),cigardel, cigardup, drpdel, drpdup, clipdel, clipdup);
		vector<vector<Signature>> SignatureDelClusters;
		vector<vector<Signature>> SignatureDupClusters;
		sortAndDeDup(ContigTypeSignatures[0]);
		clustering(ContigTypeSignatures[0],SignatureDelClusters,AllStats[i]);
		// exit(0);
		sortAndDeDup(ContigTypeSignatures[1]);
		clustering(ContigTypeSignatures[1],SignatureDupClusters,AllStats[i]);
		vector<VCFRecord> Records;
		for (int j=0;j<SignatureDelClusters.size();++j)
		{
			Records.push_back(VCFRecord(Contigs[i],Ref,SignatureDelClusters[j]));
		}
		// continue;
		for (int j=0;j<SignatureDupClusters.size();++j)
		{
			Records.push_back(VCFRecord(Contigs[i],Ref,SignatureDupClusters[j]));
		}
		sort(Records.data(),Records.data()+Records.size());
		for (int j=0;j<Records.size();++j)
		{
			if (Records[j].Keep) printf("\n%s",string(Records[j]).c_str());
		}
		//vector<Variant> ContigVariants;
		//VariantsByContig.push_back(ContigVariants);
		//callVariants(Contigs[i],VariantsByContig[VariantsByContig.size()-1],ContigSignatures,Args);
	}

	//report(VariantsByContig);
	fai_destroy(Ref);
	free(Contigs);
    return 0;
}