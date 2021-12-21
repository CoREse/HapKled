#include "clustering.h"
#include <math.h>
#include "crelib/crelib.h"
#include <vector>
using namespace std;
using namespace cre;

#define CUTE_VER
int precisionLevel(const Signature &A)
{
    if (A.Tech==1 && A.Type==1) return 0;//drp sig, imprecise pricision
    if (A.Tech==0) return 1;//SMRT CIGAR and CLIP, vague precision
    return 2;//NGS CIGAR and clip, precise precision
}

int bestPrecision(const Signature &A,const Signature &B)
{
    return MAX(precisionLevel(A),precisionLevel(B));
}

int worstPrecision(const Signature &A,const Signature &B)
{
    return MIN(precisionLevel(A),precisionLevel(B));
}

float calcOverlap(float B1, float E1, float B2, float E2)
{
    if (B1<=E2 && B1>=B2)
        return MIN(E2-B1,E1-B1);
    if (E1>=B2 && E1<=E2)
        return E1-B2;
    if (B2>=B1 && B2<=E1)//1 covers 2
        return E2-B2;
    else
        return -1;
}

/*
        PDN=max((E1.End-E1.Begin),(E2.End-E2.Begin))
        SDW=0.3
        PDW=0.2
        ODW=0.5
        SampleD=1000# if E1.Sample==E2.Sample else 2
        SD=abs((E1.End-E1.Begin)-(E2.End-E2.Begin))/max((E1.End-E1.Begin),(E2.End-E2.Begin))
        PD=min(abs(E1.Begin-E2.Begin),abs(E1.End-E2.End))
        PD=min(PD,abs((E1.Begin+E1.End)/2-(E2.Begin+E2.End)/2))
        PD/=PDN
        Overlap=calcOverlap(E1.Begin,E1.End,E2.Begin,E2.End)
        OD=Overlap/(max((E1.End-E1.Begin),(E2.End-E2.Begin)))
        OD=1-OD
        return (SD*SDW+PDW*PD+ODW*OD)*SampleD
*/

float distance(const Signature &A, const Signature &B, bool Partial, float * PPD, Stats BamStats)
{
    #ifdef CUTE_VER
    *PPD=abs(A.Begin-B.Begin);
    return *PPD;
    #endif
    //if (A.SupportedSV!=B.SupportedSV) return 100;//it shouldn't happen
    float SDW, PDW, ODW;
    float PDN=900;//TODO: should be realated by tech and type, should make pd>0.5 means not the same
    int BP=bestPrecision(A,B),WP=worstPrecision(A,B);
    if (WP==2) PDN=30;
    else if (WP==0 && BP!=1) PDN=2*(BamStats.Mean+3*BamStats.SD);//should be related with the drp stats
    else if (BP==1 && WP==0) PDN=MAX(PDN,2*(BamStats.Mean+3*BamStats.SD));
    // SDW=0.33;
    // PDW=0.33;
    // ODW=0.33;
    SDW=1;
    PDW=1;
    ODW=1;
    float SD,PD,OD;
    PD=MIN(abs(A.Begin-B.Begin),abs(A.End-B.End));
    PD=MIN(PD,abs((A.Begin+A.End)/2-(B.Begin+B.End)/2));
    PD/=PDN;
    float Overlap=calcOverlap(A.Begin,A.End,B.Begin,B.End);
    if (!Partial || precisionLevel(A)==precisionLevel(B))
    {
        SD=abs(A.Length-B.Length)/MAX(A.Length,B.Length);
        OD=Overlap/(MAX((A.End-A.Begin),(B.End-B.Begin)));
        OD=1-OD;
    }
    else
    {
        const Signature * HighPre, *LowPre;
        if (precisionLevel(A)>precisionLevel(B)) {HighPre=&A;LowPre=&B;}
        else {HighPre=&B;LowPre=&A;}
        SD=abs(HighPre->Length-LowPre->Length)/(HighPre->Length);
        OD=Overlap/(HighPre->End-HighPre->Begin);
        OD=abs(1-OD);
    }
    if (PPD!=NULL) *PPD=PD;
    if (SD>1) SD=10;
    if (PD>1) PD=10;
    if (OD>1) OD=10;
    return (SD*SDW+PDW*PD+ODW*OD);//+2.0-float(BP);
}

inline int first0(short * A, int S, int B=0)
{
    for (int i=B;i<S;++i) if (A==0) return i;
}

void simpleClustering(vector<Signature> & SortedSignatures, vector<vector<Signature>> &Clusters, Stats BamStats)//like jcrd and cuteSV, SortedSignatures may have deleted ones marked by Type=-1
{
    #ifdef CUTE_VER
    float MaxDis=200;
    int MinSupport=10;
    for (int i=0;i<SortedSignatures.size();++i)
    {
        if (SortedSignatures[i].Type==-1) continue;
        // printf("%d %d %s\n",SortedSignatures[i].Begin, SortedSignatures[i].Length, SortedSignatures[i].TemplateName.c_str());
        if (Clusters.size()==0)
        {
            Clusters.push_back(vector<Signature>());
        }
        if (Clusters.back().size()==0) Clusters.back().push_back(SortedSignatures[i]);
        else
        {
            if (SortedSignatures[i].Begin-Clusters.back().back().Begin>MaxDis)
            {
                if (Clusters.back().size()<MinSupport) Clusters.pop_back();
                Clusters.push_back(vector<Signature>());
            }
            Clusters.back().push_back(SortedSignatures[i]);
        }
    }
    if (Clusters.size()>0 && Clusters.back().size()<MinSupport) Clusters.pop_back();
    #else
    float MaxDis=0.7;
    int Size=SortedSignatures.size();
    short *Clustered=(short*) calloc(sizeof(short),Size);
    int F0=0;
    while (F0<Size)
    {
        for (int i=F0;i<SortedSignatures.size();++i)
        {
            if (Clustered[i]) continue;
            if (Clusters.size()==0)
            {
                Clusters.push_back(vector<Signature>());
            }
            bool AllLink=true;
            for (int j=0;j<Clusters.back().size();++j)
            {
                float PD;
                if (distance(SortedSignatures[i],Clusters.back()[j],true,&PD, BamStats)>MaxDis)//+bestPrecision(SortedSignatures[i],Clusters.back()[j])-2>MaxDis)
                {
                    AllLink=false;
                    if (PD+bestPrecision(SortedSignatures[i],Clusters.back()[j])-2>MaxDis)
                    {
                        Clusters.push_back(vector<Signature>());
                        Clusters.back().push_back(SortedSignatures[i]);
                        break;
                    }
                }
            }
            if (AllLink)
            {
                Clusters.back().push_back(SortedSignatures[i]);
                Clustered[i]=1;
            }
        }
        F0=first0(Clustered,Size,F0);
    }
    free(Clustered);
    #endif
    // for (int i=0;i<Clusters.size();++i)
    // {
    //     printf("%d:",Clusters[i].size());
    //     for (int j=0;j<Clusters[i].size();++j)
    //     {
    //         printf("%d %d %s,",Clusters[i][j].Begin,Clusters[i][j].Length,Clusters[i][j].TemplateName.c_str());
    //     }
    // printf("\n");
    // }
    // exit(0);
}