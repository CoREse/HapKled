#ifndef KLED_CONTIG_H
#define KLED_CONTIG_H

#include <string>

struct Contig
{
	const unsigned int ID;//Contig id in reference file
	const unsigned int Size;
	std::string Name;
	Contig(const unsigned int ID,const char *Name,unsigned int Size);
};

#endif