#ifndef KLED_KLED_H
#define KLED_KLED_H

#include <vector>

struct Arguments {
	const char * ReferenceFileName=0;
	std::vector<const char *> BamFileNames;
	std::vector<const char *> CallingContigs;
	std::string SampleName="*";
	int ThreadN=8;
	int MinSVLen=30;
	int MinMappingQuality=20;
	int MinTemplateLength=500;
	int DelMinMaxMergeDis=500;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	double DelMaxMergePortion=0.2;
	int CoverageWindowSize=100;
	int InsClipTolerance=10;
	int InsMaxGapSize=50000;
	double ASSBases[2]={10,3};//Layers of base filter scores, for average supporting scores(number of templates and signatures)
	double ASSCoverageMulti[2]={0.5,0.3};
	double LSDRSs[2]={0, 60};//For Legnth Standard Deviation Ratio Scores
	double InsASSBases[2]={10,1};
	double InsASSCoverageMulti[2]={0.5,0.3};
	double InsLSDRSs[2]={10, 60};
};

#endif