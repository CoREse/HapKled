#ifndef KLED_SIGNATURE_H
#define KLED_SIGNATURE_H

#include <string>
#include <vector>

struct Segment
{
    int Begin, End;
    Segment(int Begin,int End);
};

class Signature
{
    public:
    int Type;//0: base class, for cigar based signatures, 1: DRPSignatures, 2: ClipSignatures
    int Tech;//0: SMRT, 1:NGS
    int SupportedSV;//0: DEL, 1: Dup
    std::string TemplateName;
    //unsigned long long ReadNum;//indicate the signature is from ReadNum-th read this software reads in this contig, count by primary alignment.
    int Begin, End;//0-based
    int CN;//default -1, 0 for del, >1 for dup if known
    std::vector<Segment> Segments;
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName);
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName, Segment Read1, Segment Read2);//for drp Signatures
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName, std::vector<Segment> Segments);//for Split reads signatures
    void setCN(int CN);
};

#endif