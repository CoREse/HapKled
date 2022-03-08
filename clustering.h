#ifndef KLED_CLUSTERING_H
#define KLED_CLUSTERING_H

#include "signature.h"
#include <vector>
#include "input.h"

float distance(const Signature &A, const Signature &B, bool Partial=true, float *PPD=NULL);
int precisionLevel(const Signature &A);
int bestPrecision(const Signature &A,const Signature &B);
int worstPrecision(const Signature &A,const Signature &B);
void clustering(std::vector<Signature> & SortedSignatures, std::vector<std::vector<Signature>> &Clusters, Stats BamStats);
void simpleClustering(std::vector<Signature> & SortedSignatures, std::vector<std::vector<Signature>> &Clusters, Stats BamStats);//like jcrd and cuteSV

#endif