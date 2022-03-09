#ifndef KLED_INPUT_H
#define KLED_INPUT_H

#include "contig.h"
#include "signature.h"
#include <vector>
#include "kled.h"
#include "htslib/htslib/sam.h"

Contig * getContigs(const char * ReferenceFN, int &NSeq, int WindowSize=100);

struct Stats
{
    float MedianIS,UpIS,BelowIS,Mean,SD;
};

struct Sam
{
	htsFile* SamFile;//the file of BAM/CRAM
	bam_hdr_t* Header;//header for BAM/CRAM file
	hts_idx_t* BamIndex;
	Sam();
	void close();
};

void collectSignatures(Contig &TheContig, std::vector<Signature> *ContigTypeSignatures, Arguments & Args, std::vector<Sam>& SamFiles, std::vector<Stats> AllStats, std::vector<int> AllTechs,const char * DataSource=0);

std::vector<Stats> getAllStats(const char * ReferenceFileName, const std::vector<const char *> & BamFileNames, std::vector<int> AllTechs);

std::vector<int> getAllTechs(Arguments & Args);

std::vector<Sam> initSam(Arguments & Args);
void closeSam(std::vector<Sam> &SamFiles);

#endif