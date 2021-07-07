#include "signature.h"
using namespace std;

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName){}

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, Segment Read1, Segment Read2)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName), Segments() {Segments.push_back(Read1), Segments.push_back(Read2);}

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, vector<Segment> Segments)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName), Segments(Segments) {}

void Signature::setCN(int cn)
{
    CN=cn;
}

Segment::Segment(int Begin,int End)
:Begin(Begin),End(End){}