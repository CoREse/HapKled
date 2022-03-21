#ifndef KLED_REPORT_H
#define KLED_REPORT_H

#include <string>
#include <vector>
#include "signature.h"
#include "htslib/htslib/faidx.h"
#include "contig.h"
#include <map>
#include "kled.h"

struct HeaderEntry
{
    std::string Class;//INFO,ALT,FORMAT
    std::string ID;
    std::string Type;
    std::string Description;
    std::string Number;
    std::string Source;
    std::string Version;
    HeaderEntry(std::string Class, std::string ID, std::string Description, std::string Number="", std::string Type="", std::string Source="", std::string Version="");
    operator std::string() const;
};

class VCFHeader
{
    std::string FileFormat;
    std::string FileDate;
    std::string Reference;
    std::vector<HeaderEntry> HeaderEntries;
    std::vector<std::string> SampleNames;
    std::vector<Contig> Contigs;
    public:
    VCFHeader(const char * Reference);
    void addHeaderEntry(const HeaderEntry & Entry);
    void addSample(const char * SampleName);
    void addContig(const Contig & TheContig);
    std::string genHeader();
};

class VCFRecord
{
    public:
    std::string CHROM;
    int Pos;//0-based reference Pos of the variant, for insertion is the pos after the insertion, otherwise is the 1st base of the variant
    std::string ID;
    std::string REF;
    std::string ALT;
    std::string QUAL;
    std::string FILTER;
    std::string INFO;
    std::map<std::string,std::string> Sample;

    bool Keep;//keep this record

    VCFRecord(const Contig & TheContig, faidx_t * Ref, std::vector<Signature> & SignatureCluster, double* CoverageWindows, double WholeCoverage, Arguments &Args, double* CoverageWindowsSums=NULL, double* CheckPoints=NULL, int CheckPointInterval=0);
    operator std::string() const;
    bool operator<(const VCFRecord& Other) const;
};

void addKledEntries(VCFHeader & Header);

double getAverageCoverage(int Begin, int End, double * CoverageWindows, Arguments & Args, double* CoverageWindowsSums=NULL, double* CheckPoints=NULL, int CheckPointInterval=0);

#endif