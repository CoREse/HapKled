#include "signature.h"
#include <string.h>
using namespace std;

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName), Length(End-Begin){}

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, Segment Read1, Segment Read2, int ALength)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName), Segments(), Length(ALength) {Segments.push_back(Read1), Segments.push_back(Read2);}

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, vector<Segment> Segments)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName), Segments(Segments), Length(End-Begin) {}

void Signature::setCN(int cn)
{
    CN=cn;
}
bool Signature::operator<(const Signature &Other) const
{
    char AL[100], BL[100];
    sprintf(AL,"%d",this->Length);
    sprintf(BL,"%d",Other.Length);
    int cmp=strcmp(AL,BL);
    if (this->Begin<Other.Begin) return true;
    else if (this->Begin>Other.Begin) return false;
    // else if (this->Length<Other.Length) return true;
    // else if (this->Length>Other.Length) return false;
    else if (cmp<0) return true;
    else if (cmp>0) return false;
    else return this->TemplateName<Other.TemplateName;
}

bool Signature::operator==(const Signature &Other) const
{
    return Tech==Other.Tech && SupportedSV==Other.SupportedSV && TemplateName==Other.TemplateName && Begin==Other.Begin && End==Other.End && Segments==Other.Segments;
}

Segment::Segment(int Begin,int End)
:Begin(Begin),End(End){}

bool Segment::operator==(const Segment & Other) const
{
    return Begin==Other.Begin && End==Other.End;
}