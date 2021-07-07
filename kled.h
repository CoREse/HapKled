#include <vector>

struct Arguments {
	const char * ReferenceFileName=0;
	std::vector<const char *> BamFileNames;
    int MinSVLen=30;
};
extern Arguments Args;
