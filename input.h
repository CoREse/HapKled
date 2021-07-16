#ifndef KLED_INPUT_H
#define KLED_INPUT_H

#include "contig.h"
#include "signature.h"
#include <vector>
#include "kled.h"

Contig * getContigs(const char * ReferenceFN, int &NSeq, int WindowSize=100);

struct Stats
{
    float MedianIS,UpIS,BelowIS,Mean,SD;
};

void collectSignatures(Contig &TheContig, std::vector<Signature> *ContigTypeSignatures, Arguments & Args, std::vector<Stats> AllStats, std::vector<int> AllTechs,const char * DataSource="samtools");

std::vector<Stats> getAllStats(const char * ReferenceFileName, const std::vector<const char *> & BamFileNames, std::vector<int> AllTechs);

std::vector<int> getAllTechs(Arguments & Args);

#endif