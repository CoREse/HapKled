#include "report.h"
#include <assert.h>
#include "clustering.h"
#include <tuple>
#include <algorithm>
#include <math.h>
#include "crelib/crelib.h"
#include "time.h"
#include <set>
#include <numeric>
#include "defines.h"

using namespace std;

float stdStat(Signature* Begin,Signature* End)
{
    float Mean=0;
    int k=0;
    for (Signature * i=Begin;i!=End;++i)
    {
        Mean+=i->Begin;
        ++k;
    }
    Mean/=float(k);
    float std=0;
    for (Signature * i=Begin;i!=End;++i)
    {
        std+=pow(i->Begin-Mean,2);
    }
    std/=float(k);
    return pow(std,0.5);
}

tuple<int,int> consensusFixed(vector<Signature> &Signatures, int K, float stat(Signature*,Signature*))//计算标准差最小的连续K个的pos和length的均值
{
    if (Signatures.size()==1) return tuple<int,int>(Signatures[0].Begin,(Signatures[0].End-Signatures[0].Begin));
    float MinStd=0;
    int MinI=-1;
    for (int i=0;i<Signatures.size()-K;++i)
    {
        float Std=stat(Signatures.data()+i,Signatures.data()+i+K);
        if (MinI==-1 || MinStd<Std)
        {
            MinI=i;
            MinStd=Std;
        }
    }
    float Pos=0,Length=0;
    for (int i=MinI;i<MinI+K;++i)
    {
        Pos+=Signatures[i].Begin;
        Length+=Signatures[i].Length;
    }
    return tuple<int,int>(int(Pos/float(K)),int(Length/float(K)));
}

const char * getSVType(const vector<Signature> &Signatures)
{
    return Signatures[0].SupportedSV==0?"DEL":"DUP";
}

tuple<int,int> getConsensusSite(vector<Signature> &Signatures)
{
    if (Signatures.size()==0) return tuple<int,int>(-1,-1);
    sort(Signatures.data(),Signatures.data()+Signatures.size());
    return consensusFixed(Signatures,Signatures.size()/2,stdStat);
}

tuple<int,int> analyzeSignatureCluster(vector<Signature> &SignatureCluster)
{
    float L2Weight=0.75;
    int Pos,Length;
    vector<Signature> PSigs[3];
    for (int i=0;i<SignatureCluster.size();++i)
    {
        PSigs[precisionLevel(SignatureCluster[i])].push_back(SignatureCluster[i]);
    }
    tuple<int,int> Sites[3];
    for (int i=0;i<3;++i) Sites[i]=getConsensusSite(PSigs[i]);
    if (get<0>(Sites[2])!=-1)//has pricise level signatures
    {
        Pos=get<0>(Sites[2]);
        Length=get<1>(Sites[2]);
    }
    else
    {
        if (get<0>(Sites[1])==-1)//doesn't have vague level signatures
        {
            Pos=get<0>(Sites[0]);
            Length=get<1>(Sites[0]);
        }
        else if (get<0>(Sites[0])==-1)//doesn't have imprecise level signatures
        {
            Pos=get<0>(Sites[1]);
            Length=get<1>(Sites[1]);
        }
        else
        {
            Pos=get<0>(Sites[0])*(1.0-L2Weight)+get<0>(Sites[1])*(L2Weight);
            Length=get<1>(Sites[0])*(1.0-L2Weight)+get<1>(Sites[1])*(L2Weight);
        }
    }
    return tuple<int,int>(Pos,Length);
}

template <typename IterType> double calcSD(IterType Begin, IterType End)//Not estimate, just calc
{
    int Size=0;
    double Mean=0;
    for (IterType it=Begin;it!=End;++it)
    {
        Mean+=*it;
        ++Size;
    }
    Mean/=(double)Size;
    double SD=0;
    for (IterType it=Begin;it!=End;++it)
    {
        SD+=pow((*it)-Mean,2.0);
    }
    SD/=(double)Size;
    return pow(SD,0.5);
}

template <typename ValueType,typename IterType> class MemberIter
{
    protected:
    IterType BaseIter;
    public:
    MemberIter(IterType A):BaseIter(A){}
    MemberIter& operator++()
    {
        ++BaseIter;
        return *this;
    }
    bool operator!=(const MemberIter & Other)
    {
        return BaseIter!=Other.BaseIter;
    }
    bool operator!=(const IterType & Other)
    {
        return BaseIter!=Other;
    }
    // virtual ValueType operator* ();//Can't use virtual in template!
};

template <typename ValueType,typename IterType> class BeginIter : public MemberIter<ValueType,IterType>
{
    public:
    // using MemberIter<ValueType,IterType>::MemberIter<ValueType,IterType>;
    BeginIter(IterType A):MemberIter<ValueType,IterType>(A){}
    ValueType operator* ()
    {
        return this->BaseIter->Begin;
    }
};

template <typename ValueType,typename IterType> class EndIter : public MemberIter<ValueType,IterType>
{
    public:
    EndIter(IterType A):MemberIter<ValueType,IterType>(A){}
    ValueType operator* ()
    {
        return this->BaseIter->End;
    }
};
template <typename ValueType,typename IterType> class LengthIter : public MemberIter<ValueType,IterType>
{
    public:
    LengthIter(IterType A):MemberIter<ValueType,IterType>(A){}
    ValueType operator* ()
    {
        return this->BaseIter->Length;
    }
};

vector<double> scoring(vector<Signature> &SignatureCluster,int SVLen)
{
    SVLen=abs(SVLen);
    vector<double> Scores;//supported reads, supported signatures, (begin sd+end sd)/2, length sd, number of signature source
    int SS=0;
    set<string> SupportTemps;;
    for (int i =0;i<SignatureCluster.size();++i)
    {
        ++SS;
        SupportTemps.insert(SignatureCluster[i].TemplateName);
    }
    int ST=SupportTemps.size();
    double STS=ST>50?100.0:double(ST)/50.0*100.0;
    double SSS=SS>50?100.0:double(SS)/50.0*100.0;
    double ASS=(STS+SSS)/2.0;
    Scores.push_back(ASS);
    // Scores.push_back(STS);
    // Scores.push_back(SSS);
    double BESD=calcSD(BeginIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),BeginIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    BESD+=calcSD(EndIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),EndIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    BESD/=2.0;
    double LengthSD=calcSD(LengthIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),LengthIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    double BESDS=BESD>50?100.0:double(BESD)/50.0*100.0;BESDS=100.0-BESDS;
    double LengthSDS=LengthSD>50?100.0:double(LengthSD)/50.0*100.0;LengthSDS=100.0-LengthSDS;
    Scores.push_back(BESDS);
    double BESDRatio=BESD/(double)SVLen;
    double BESDRatioS=BESDRatio>0.5?100:BESDRatio/0.5*100.0;BESDRatioS=100.0-BESDRatioS;
    Scores.push_back(BESDRatioS);
    Scores.push_back(LengthSDS);
    double LengthSDRatio=LengthSD/(double)SVLen;
    double LengthSDRatioS=LengthSDRatio>0.5?100:LengthSDRatio/0.5*100.0;LengthSDRatioS=100.0-LengthSDRatioS;
    Scores.push_back(LengthSDRatioS);
    set<int> Sources;
    for (int i =0;i<SignatureCluster.size();++i)
    {
        if (SignatureCluster[i].Type!=-1) Sources.insert(SignatureCluster[i].Type);
    }
    Scores.push_back(double(Sources.size())/3.0*100.0);
    return Scores;
}

double getAmbientCoverage(int Begin, int End, double * CoverageWindows, Arguments & Args)
{
    if (End<=Begin) return -1;
    int WBegin=Begin/Args.CoverageWindowSize;
    int WEnd=End/Args.CoverageWindowSize+1;
    double Cov=0;
    for (int i=WBegin;i<WEnd;++i) 
    {
        Cov*=((double)(i-WBegin))/(i-WBegin+1);//in case too big
        Cov+=CoverageWindows[i]*1/(i-WBegin+1);
    }
    return Cov;
}

bool keepCluster(vector<Signature> &SignatureCluster, int & SS, int &ST)
{
    SS=0;
    set<string> SupportTemps;;
    for (int i =0;i<SignatureCluster.size();++i)
    {
        ++SS;
        SupportTemps.insert(SignatureCluster[i].TemplateName);
    }
    ST=SupportTemps.size();
    if (ST>=2) return true;
    return false;
    double LengthSD=calcSD(LengthIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),LengthIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    // if (ST>=12 && LengthSD<20) return true;
    if (ST>=10 && LengthSD<10) return true;
    if (ST>=5 && LengthSD<5 ) return true;
    set<int> Sources;
    for (int i =0;i<SignatureCluster.size();++i)
    {
        if (SignatureCluster[i].Type!=-1) Sources.insert(SignatureCluster[i].Type);
    }
    // if (ST>=3 && LengthSD<5 && Sources.size()>1) return true;
    return false;
}

int VN=0;

VCFRecord::VCFRecord(const Contig & TheContig, faidx_t * Ref,vector<Signature> & SignatureCluster, double* CoverageWindows, double WholeCoverage, Arguments& Args)
{
    // printf("%d:",SignatureCluster.size());
    // for (int j=0;j<SignatureCluster.size();++j)
    // {
    //     printf("%d %d %s,",SignatureCluster[j].Begin,SignatureCluster[j].Length,SignatureCluster[j].TemplateName.c_str());
    // }
    // printf("\n");
    // return;
    assert(SignatureCluster.size()>=0);
    int SS,ST;
    // printf("%d, ",SignatureCluster.size());
    // for (int i=0;i<SignatureCluster.size();++i)
    // {
    //     // printf("[%d, %d, '%s'] ",SignatureCluster[i].Begin,SignatureCluster[i].Length,SignatureCluster[i].TemplateName.c_str());
    //     printf("%d ", SignatureCluster[i].Length);
    // }
    // printf("\n");
    // return;
    if (keepCluster(SignatureCluster,SS,ST)) Keep=true;
    else {Keep=false; return;}
    string SVType=getSVType(SignatureCluster);
    tuple<int,int> Site=analyzeSignatureCluster(SignatureCluster);
    int SVLen=get<1>(Site);
    Pos=get<0>(Site);//0-bsed now, after ref and alt then transform to 1-based, but should be the base before variantion. End should be the last base, but also should be transform to 1-based. So they don't change.
    //VCF version 4.2 says if alt is <ID>, the pos is the base preceding the polymorphism. No mention of the "base 1" rule.
    #ifdef CUTE_VER
    unsigned long long llPos=0;
    unsigned long long llSVLen=0;
    for (int i=0;i<SignatureCluster.size();++i)
    {
        llPos+=SignatureCluster[i].Begin;
        llSVLen+=SignatureCluster[i].Length;
    }
    Pos=llPos/SignatureCluster.size();
    SVLen=llSVLen/SignatureCluster.size();
    #endif
    unsigned long long llPos=0;
    unsigned long long llSVLen=0;
    for (int i=0;i<SignatureCluster.size();++i)
    {
        llPos+=SignatureCluster[i].Begin;
        llSVLen+=SignatureCluster[i].Length;
    }
    Pos=llPos/SignatureCluster.size();
    SVLen=llSVLen/SignatureCluster.size();

    
    vector<double> Scores=scoring(SignatureCluster,SVLen);
    double ScoreWeights[6]={0.29279039862777806, 0.015320183836380931, 0.14398052205008294, 0.17979354517797344, 0.2617766118686123, 0.10633873843917234};
    double Score=0;for (int i=0;i<Scores.size();++i) Score+=ScoreWeights[i]*Scores[i];
    double AmbientCoverage=getAmbientCoverage(Pos,Pos+abs(SVLen),CoverageWindows,Args);
    // if (Score<60) {Keep=false;return;}
    // if (Scores[0]>=3 && Scores[4]>=55 && AmbientCoverage<WholeCoverage*0.6) Keep=true;
    if (Scores[0]>10+WholeCoverage*0.5) {Keep=true;}
    else
    {
        if (Scores[0]<3+WholeCoverage*0.3) {Keep=false;return;}
        // if (Scores[1]<95) {Keep=false;return;}
        if (Scores[4]<60) {Keep=false;return;}
        // if (abs(SVLen)>10000)
        // {
        //     double BeforeCovergae=-1;
        //     if (Pos>10000) BeforeCovergae=getAmbientCoverage(Pos-10000,Pos,CoverageWindows,Args);
        //     double AfterCoverage=-1;
        //     if (Pos+abs(SVLen)<TheContig.Size-10000) AfterCoverage=getAmbientCoverage(Pos+abs(SVLen),Pos+abs(SVLen)+10000,CoverageWindows,Args);
        //     double Ambient=0;
        //     double Valid=0;
        //     if (BeforeCovergae!=-1) {Valid+=1;Ambient+=BeforeCovergae;}
        //     if (AfterCoverage!=-1) {Valid+=1;Ambient+=AfterCoverage;}
        //     if (Valid!=0) Ambient/=Valid;
        //     if (Ambient!=0)
        //     {
        //         double Center=getAmbientCoverage(Pos+abs(SVLen)/4.0,Pos+abs(SVLen)*0.75,CoverageWindows,Args);
        //         if (Center>Ambient*0.6) {Keep=false;return;}
        //     }
        // }
        // if (abs(SVLen)>5000 && AmbientCoverage>WholeCoverage*1.0) {Keep=false;return;}
        // else if (Scores[0]>=10 && (Scores[2]+Scores[3]+Scores[4])>=240) Keep=true;
        // else if (Scores[0]>=6 && ((Scores[2]+Scores[3])>=190 || (Scores[2]+Scores[4])>=190 || (Scores[3]+Scores[4])>=190 )) Keep=true;
        // else {Keep=false;return;}
    }
    int TLen;
    int End;
    bool OutTag=true;
    if (SVLen>10000) OutTag=true;
    if (OutTag)
    {
        if (Pos==0)
        {
            Pos=SVLen;
            End=Pos-1;
            char * TSeq=faidx_fetch_seq(Ref,TheContig.Name.c_str(),Pos,Pos,&TLen);
            REF=TSeq[0];
            free(TSeq);
            ALT="<"+SVType+">"+REF;
        }
        else
        {
            --Pos;//base before variant
            End=Pos+SVLen;
            char * TSeq=faidx_fetch_seq(Ref,TheContig.Name.c_str(),Pos,Pos,&TLen);
            REF=TSeq[0];
            free(TSeq);
            ALT="<"+SVType+">";
        }
    }
    else
    {
        if (Pos==0)
        {
            Pos=SVLen;
            End=Pos-1;
            char * TSeq=faidx_fetch_seq(Ref,TheContig.Name.c_str(),0,Pos,&TLen);
            REF=string(TSeq);
            free(TSeq);
            if (SVType=="DEL")
                ALT=REF[Pos];
            else
                ALT="<"+SVType+">"+REF[Pos];
        }
        else
        {
            --Pos;//base before variant
            End=Pos+SVLen;
            char * TSeq=faidx_fetch_seq(Ref,TheContig.Name.c_str(),Pos-1,End,&TLen);
            REF=TSeq;
            free(TSeq);
            if (SVType=="DEL")
                ALT=REF[0];
            else
                ALT="<"+SVType+">";
        }
    }
    ++Pos;++End;//trans to 1-based
    INFO="PRECISE;SVTYPE="+SVType+";END="+to_string(End)+";SVLEN="+to_string(SVType=="DEL"?-SVLen:SVLen)+";SS="+to_string(SS)+";ST="+to_string(ST);
    INFO+=";SCORES="+to_string(int(Scores[0]));
    for (int i=1;i<Scores.size();++i) INFO+=","+to_string(int(Scores[i]));
    //extern int VN;
    //ID="kled."+SVType+"."+to_string(VN);
    //++VN;
    ID=".";
    QUAL=".";
    FILTER="PASS";
    CHROM=TheContig.Name;
    Sample["GT"]="./.";
}

VCFRecord::operator std::string() const
{
    string Output=CHROM+"\t"+to_string(Pos)+"\t"+ID+"\t"+REF+"\t"+ALT+"\t"+QUAL+"\t"+FILTER+"\t"+INFO+"\t";
    bool FirstFormat=true;
    string SampleString=("");
    for (auto i=Sample.begin();i!=Sample.end();++i)
    {
        if (!FirstFormat) Output+=":";
        Output+=i->first;
        SampleString+=i->second;
        FirstFormat=false;
    }
    Output+="\t"+SampleString;
    return Output;
}

bool VCFRecord::operator<(const VCFRecord& Other) const
{
    return Pos<Other.Pos;
}

HeaderEntry::HeaderEntry(string AClass, string AID, string ADescription, string ANumber, string AType, string ASource, string AVersion)
:Class(AClass),ID(AID),Description(ADescription),Type(AType),Number(ANumber),Source(ASource),Version(AVersion)
{
}
HeaderEntry::operator string() const
{
    if (Class=="INFO" || Class=="FORMAT")
        return "##"+Class+"=<ID="+ID+",Number="+Number+",Type="+Type+",Description=\""+Description+"\""+(Source!=""?(",Source=\""+Source+"\""):"")+(Version!=""?",Version=\""+Version+"\"":"")+">";
    else
        return "##"+Class+"=<ID="+ID+",Description=\""+Description+"\">";
}

inline string getDate()
{
    time_t t = time(0);   // get time now
    tm* now = localtime(&t);
    return to_string(now->tm_year + 1900)+(now->tm_mon<9?"0":"")+to_string((now->tm_mon + 1))+to_string(now->tm_mday);
}

VCFHeader::VCFHeader(const char * AReference)
{
    FileFormat="##fileformat=VCFv4.2";
    FileDate="##fileDate="+getDate();
    Reference=string("##reference=")+string(AReference);
}
void VCFHeader::addHeaderEntry(const HeaderEntry & Entry)
{
    HeaderEntries.push_back(Entry);
}
void VCFHeader::addSample(const char * SampleName)
{
    SampleNames.push_back(SampleName);
}

void VCFHeader::addContig(const Contig & TheContig)
{
    Contigs.push_back(TheContig);
}

std::string VCFHeader::genHeader()
{
    string Output;
    Output+=FileFormat+"\n"+FileDate+"\n"+Reference;
    for (int i=0;i<Contigs.size();++i)
    {
        Output+=string("\n")+"##contig=<ID="+Contigs[i].Name+",length="+to_string(Contigs[i].Size)+">";
    }
    for (int i=0;i<HeaderEntries.size();++i)
    {
        Output+=string("\n")+string(HeaderEntries[i]);
    }
    Output+="\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int i=0;i<SampleNames.size();++i) Output+="\t"+SampleNames[i];
    return Output;
}

void addKledEntries(VCFHeader & Header)
{
    Header.addHeaderEntry(HeaderEntry("INFO","SVTYPE","Type of structural variant","1","String"));
    Header.addHeaderEntry(HeaderEntry("INFO","END","End position of the variant described in this record","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","SVLEN","Difference in length between REF and ALT alleles","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","PRECISE","Precise structural variant","0","Flag"));
    Header.addHeaderEntry(HeaderEntry("INFO","IMPRECISE","Imprecise structural variant","0","Flag"));
    Header.addHeaderEntry(HeaderEntry("INFO","SS","Number of supported signatures","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","ST","Number of supported templates","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("ALT","DEL","Deletion"));
    Header.addHeaderEntry(HeaderEntry("ALT","DUP","Duplication"));
    Header.addHeaderEntry(HeaderEntry("FORMAT","GT","Genotype","1","String"));
    Header.addHeaderEntry(HeaderEntry("INFO","SCORES","Scores, supported reads, supported signatures, (begin sd+end sd)/2, length sd, number of signature source","6","Integer"));
}