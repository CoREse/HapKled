#include <string>
#include <vector>
#include "signature.h"

struct AlignmentSigs
{
	unsigned long long AlignmentID;
	std::string TemplateName;
	std::vector<std::vector<Signature>> TypeSignatures;
	int BeginMost, EndMost;
	std::vector<int> TypeBeginMost, TypeEndMost;
	AlignmentSigs(unsigned long long AlignmentID, const char *TempName="");
	int getBeginMost();
	int getEndMost();
};

struct MergingMutex
{
	pthread_mutex_t m_AddingSigs[2];
	MergingMutex()
	{
		pthread_mutex_init(&m_AddingSigs[0],NULL);
		pthread_mutex_init(&m_AddingSigs[1],NULL);
	}
};

struct OmniBMergeArgs
{
	std::vector<std::vector<Signature>> * pTypeSignatures;
	int Index;
	std::vector<AlignmentSigs> * pAlignmentsSigs;
	const int ** TypeSigIndexes;
	const int ** TypeMaxEnds;
	MergingMutex *mut;
    Arguments *pArgs;
};

void *omniBHandler(void * Args);
void divideASs(std::vector<AlignmentSigs> &ASs);

struct SimpleMergeArgs
{
	std::vector<std::vector<Signature>> * pTypeSignatures;
    int Index;
    std::vector<AlignmentSigs> * pAlignmentsSigs;
    MergingMutex *mut;
    Arguments *pArgs;
	bool Regional;
};

void *simpleMergeHandler(void * Args);