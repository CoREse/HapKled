#ifndef KLED_VARIANT_H
#define KLED_VARIANT_H

#include <string>

class Variant
{
    public:
    int Type;
    std::string ContigName;
    int Pos;
    int End;
    int SVLen;
};

#endif