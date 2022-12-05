#pragma once
#ifndef KLED_KLED_H
#define KLED_KLED_H

#include <vector>

const int NumberOfSVTypes=4;//Default is static so is fine.

struct Arguments {
	const char * Version="1.0";
	bool ShowVersion=false;
	bool ShowHelp=false;
	std::string CommandLine;
	std::string RunHash;//Storing std::hash<std::string>(Version+Arguments) in hex.
	int TestN=0;
	const char * ReferenceFileName=0;
	std::vector<const char *> BamFileNames;
	std::vector<const char *> CallingContigs;
	std::string SampleName="*";
	int ThreadN=8;
	bool NoFilter=false;
	bool Filter2ST=false;
	int MinSVLen=30;
	int MinMappingQuality=20;
	int MinTemplateLength=500;
	int DelMinMaxMergeDis=500;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	int DelMaxMaxMergeDis=50000;
	int InsMinMaxMergeDis=500;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	double DelMaxMergePortion=0.2;
	int CoverageWindowSize=100;
	int InsClipTolerance=10;
	int InsMaxGapSize=50000;
	int MaxClusterSize=1000;//will limit (sampling) Cluster size to [MaxClusterSize, 2*MaxClusterSize)
	int ClusteringMaxMergeRange=20000;
	int ClusteringBatchSize=10000;
    float BrotherhoodTypeRatios[NumberOfSVTypes]={0.5,1.0,0.1,0.5};//For SV Types
    int BrotherhoodTypeForceBrothers[NumberOfSVTypes]={10,50,500,100};
    float BrotherhoodTypeLengthRatios[NumberOfSVTypes]={0.5,0.3,0.1,0.5};//For SV Types
    int BrotherhoodTypeLengthMinEndurance[NumberOfSVTypes]={10,0,500,100};
    int BrotherhoodNearRanges[NumberOfSVTypes]={-1,-1,200,200};//If -1, use brotherhoodcluster1
    int BrotherhoodTypeForceBrothers2[NumberOfSVTypes]={10,50,500,100};
    float BrotherhoodTypeLengthRatios2[NumberOfSVTypes]={0.5,0.3,0.1,0.5};//For SV Types
    float BrotherhoodCCSTypeRatios[NumberOfSVTypes]={2.0,2.5,2.0,2.0};//For SV Types
    int BrotherhoodCCSTypeForceBrothers[NumberOfSVTypes]={100,200,100,100};
    float BrotherhoodCCSTypeLengthRatios[NumberOfSVTypes]={0.5,0.5,0.1,0.5};//For SV Types
    int BrotherhoodCCSTypeLengthMinEndurance[NumberOfSVTypes]={10,50,500,100};
	bool AllCLR=false;
	bool AllCCS=false;
	double ASSBases[NumberOfSVTypes][2]=//Layers of base filter of the addition of supporting segmentations and templates
	{{1,1}//DEL
	,{3,4}//INS
	,{0,2}//DUP
	,{0,2}};//INV, only ST no SS
	double ASSCoverageMulti[NumberOfSVTypes][2]={{0.2,0.1},
	{0,0},
	{0.4,0.4},
	{0.2,0.3}};
	double LSDRSs[NumberOfSVTypes][2]={{0, 50},//For Legnth Standard Deviation Ratio Scores(100-ratio*100)
	{40,0},
	{95,90},
	{94,70}};
	double CLRASSBases[NumberOfSVTypes][2]=//Layers of base filter of the addition of supporting segmentations and templates
	{{1.5,0.5}//DEL
	,{2,2}//INS
	,{0,2}//DUP
	,{1,1}};//INV, only ST no SS
	double CLRASSCoverageMulti[NumberOfSVTypes][2]={{0.1,0.1},
	{0.18,0.08},
	{0.4,0.4},
	{0.5,0.5}};
	double CLRLSDRSs[NumberOfSVTypes][2]={{30, 95},//For Legnth Standard Deviation Ratio Scores(100-ratio*100)
	{55,85},
	{95,90},
	{55,90}};
	double CCSASSBases[NumberOfSVTypes][2]={{10,3},
	{5,1},
	{10,3},//Need further polish
	{1,1}};//Need further polish
	double CCSASSCoverageMulti[NumberOfSVTypes][2]={{0.5,0.1},
	{0.3,0.2},
	{0.5,0.1},
	{0.5,0.5}};
	double CCSLSDRSs[NumberOfSVTypes][2]={{0, 80},
	{55,80},
	{0,80},
	{55,90}};
    double PreciseStandard=3;
    int MinimumPreciseTemplates=5;

	double TotalCoverage=0;

	std::string CustomFilterParas[NumberOfSVTypes]={"","","",""};
	std::string CustomClusterParas[NumberOfSVTypes]={"","","",""};

	int CallByContig=false;
};

#endif
