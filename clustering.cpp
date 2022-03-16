#include "clustering.h"
#include <math.h>
#include "crelib/crelib.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "defines.h"
#include <list>
using namespace std;
using namespace cre;

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

bool SigLengthLess(const Signature & a, const Signature &b)
{
    if (a.Length<b.Length) return true;
    if (a.Length>b.Length) return false;
    if (a.Begin<b.Begin) return true;
    if (a.Begin>b.Begin) return false;
    return a.TemplateName<b.TemplateName;
}

void keepLongestPerRead(vector<Signature> & SignatureCluster)
{
    unordered_map<string, Signature> TempsSig;
    for (int i=0;i<SignatureCluster.size();++i)
    {
        if(TempsSig.count(SignatureCluster[i].TemplateName)==0)
        {
            TempsSig[SignatureCluster[i].TemplateName]=SignatureCluster[i];
        }
        else if (TempsSig[SignatureCluster[i].TemplateName].Length<SignatureCluster[i].Length)
        {
            TempsSig[SignatureCluster[i].TemplateName]=SignatureCluster[i];
        }
    }
    SignatureCluster.clear();
    for (auto it=TempsSig.begin();it!=TempsSig.end();++it)
    {
        SignatureCluster.push_back(it->second);
    }
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
    //add dedup and length segment here
    vector<vector<Signature>> OldCs=Clusters;
    Clusters.clear();
    for (int i=0;i< OldCs.size();++i)
    {
        keepLongestPerRead(OldCs[i]);
        sort(OldCs[i].begin(),OldCs[i].end(),SigLengthLess);
        if (OldCs[i].size()<10) continue;
        double MeanLength=0;
        for (int j=0;j<OldCs[i].size();++j)
        {
            MeanLength+=OldCs[i][j].Length;
        }
        MeanLength/=double(OldCs[i].size());
        double Discrete=0.5*MeanLength;
        Clusters.push_back(vector<Signature>());
        Clusters.back().push_back(OldCs[i][0]);
        int LastLen=OldCs[i][0].Length;
        for (int j=1;j<OldCs[i].size();++j)
        {
            if (OldCs[i][j].Length-LastLen>Discrete)
            {
                Clusters.push_back(vector<Signature>());
            }
            Clusters.back().push_back(OldCs[i][j]);
            LastLen=OldCs[i][j].Length;
        }
    }
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
    //         printf("%d,",Clusters[i][j].Length);
    //     }
    // printf("\n");
    // }
    // exit(0);
}

bool isBrother(const Signature &A, const Signature &B, float Ratio=0.1, int ForceBrother=5)
{
    if (abs(A.Begin-B.Begin)<=ForceBrother && abs(A.End-B.End)<=ForceBrother) return true;
    int MinLength=min(A.Length,B.Length);
    if (abs(A.Begin-B.Begin)<=MinLength*Ratio && abs(A.End-B.End)<=MinLength*Ratio && abs(A.Length-B.Length)<=MinLength*Ratio) return true;
    return false;
}

class Brotherhood
{
    public:
    vector<Signature> Cluster;
    bool CCS;
    unsigned int MinBegin, MaxEnd;
    inline static float TypeRatios[3]={0.5,0.5,0.5};//For SV Types
    inline static int TypeForceBrothers[3]={10,50,10};
    inline static float CCSTypeRatios[3]={1.5,1.5,0.5};//For SV Types
    inline static int CCSTypeForceBrothers[3]={100,100,10};
    Brotherhood(bool IsCCS=false):Cluster(),MinBegin(0),MaxEnd(0),CCS(IsCCS) {};
    Brotherhood(Signature &S, bool IsCCS=false):Cluster(),MinBegin(0),MaxEnd(0),CCS(IsCCS) {setCluster(S);};
    void setCluster(Signature & S)
    {
        Cluster.clear();
        Cluster.push_back(S);
        MinBegin=S.Begin;
        MaxEnd=S.End;
    }
    bool canMerge(const Brotherhood & Other) const
    {
        int ForceBrother=Brotherhood::TypeForceBrothers[Cluster[0].SupportedSV];
        float Ratio=Brotherhood::TypeRatios[Cluster[0].SupportedSV];
        if (CCS)
        {
            ForceBrother=Brotherhood::CCSTypeForceBrothers[Cluster[0].SupportedSV];
            Ratio=Brotherhood::CCSTypeRatios[Cluster[0].SupportedSV];
        }
            // fprintf(stderr,"%d %d %d %d %d %d %d\n",Other.MinBegin-MaxEnd, Other.MinBegin, MaxEnd, MinBegin-Other.MaxEnd,MinBegin, Other.MaxEnd,Brotherhood::ForceBrother);
        if (((Other.MinBegin>MaxEnd+ForceBrother) || ((MinBegin>Other.MaxEnd+ForceBrother)))) return false;
            // fprintf(stderr,"yes");
        for (const Signature & A:Cluster)
        {
            for (const Signature & B:Other.Cluster)
            {
                if (isBrother(A,B,Ratio,ForceBrother))
                {
                    return true;
                }
            }
        }
        return false;
    }
    bool merge(Brotherhood & Other)
    {
        if (canMerge(Other))
        {
            for (Signature & B:Other.Cluster)
            {
                Cluster.push_back(B);
            }
            MinBegin=min(MinBegin,Other.MinBegin);
            MaxEnd=max(MaxEnd,Other.MaxEnd);
            return true;
        }
        return false;
    }
};

void brotherClustering(vector<Signature> & SortedSignatures, vector<vector<Signature>> &Clusters, Stats BamStats, Arguments &Args,int MaxMergeRange=20000)
{
    list<Brotherhood> Brotherhoods;
    for (Signature & S:SortedSignatures) Brotherhoods.push_back(Brotherhood(S,Args.AllCCS));
    while (1)
    {
        Brotherhoods.sort([](Brotherhood & a, Brotherhood &b)-> bool {return a.MinBegin<b.MinBegin;});
        bool Next=false;
        for (list<Brotherhood>::iterator Ai=Brotherhoods.begin();Ai!=Brotherhoods.end();++Ai)
        {
            for (list<Brotherhood>::iterator Bi=next(Ai);Bi!=Brotherhoods.end();++Bi)
            {
                if (Ai->MinBegin+MaxMergeRange<Bi->MinBegin) break;
                if (Ai->merge(*Bi))
                {
                    Brotherhoods.erase(Bi);
                    Next=true;
                    break;
                }
            }
        }
        if (!Next) break;
    }
    for (Brotherhood& B :Brotherhoods)
    {
        Clusters.push_back(B.Cluster);
    }
}

void clustering(vector<Signature> & SortedSignatures, vector<vector<Signature>> &Clusters, Stats BamStats, Arguments& Args)
{
    fprintf(stderr,"%lu %lu:\n",SortedSignatures.size(), Clusters.size());
    brotherClustering(SortedSignatures,Clusters,BamStats,Args);
    // simpleClustering(SortedSignatures,Clusters,BamStats);
    fprintf(stderr,"%lu %lu:\n",SortedSignatures.size(), Clusters.size());
    // for (int i=0;i<Clusters.size();++i) fprintf(stderr, "%d\n",Clusters[i].size());
}