#include "signature.h"
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
bool Signature::operator<(const Signature &Other)
{
    return this->Begin<Other.Begin;
}

Segment::Segment(int Begin,int End)
:Begin(Begin),End(End){}