#ifndef KLED_INPUT_H
#define KLED_INPUT_H

#include "contig.h"
#include "signature.h"
#include <vector>

Contig * getContigs(const char * ReferenceFN, int &NSeq, int WindowSize=100);

struct Stats
{
    float MedianIS,UpIS,BelowIS,Mean,SD;
};

void collectSignatures(Contig &TheContig, std::vector<Signature> &ContigSignatures, const char * ReferenceFileName, const std::vector<const char *> & BamFileNames, std::vector<Stats> AllStats, const char * DataSource="samtools");

std::vector<Stats> getAllStats(const char * ReferenceFileName, const std::vector<const char *> & BamFileNames);

#endif