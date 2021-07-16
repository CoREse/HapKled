#include "report.h"
#include <assert.h>
#include "clustering.h"
#include <tuple>
#include <algorithm>
#include <math.h>
#include "crelib/crelib.h"
#include "time.h"
#include <set>

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

tuple<int,int> consensusFixed(vector<Signature> &Signatures, int K, float stat(Signature*,Signature*))
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

int VN=0;

bool keepCluster(vector<Signature> SignatureCluster, int & SS, int &ST)
{
    SS=0;
    set<string> SupportTemps;;
    for (int i =0;i<SignatureCluster.size();++i)
    {
        ++SS;
        SupportTemps.insert(SignatureCluster[i].TemplateName);
    }
    ST=SupportTemps.size();
    if (ST>=5) return true;
    return false;
}

VCFRecord::VCFRecord(const Contig & TheContig, faidx_t * Ref,vector<Signature> & SignatureCluster)
{
    assert(SignatureCluster.size()>=0);
    int SS,ST;
    if (keepCluster(SignatureCluster,SS,ST)) Keep=true;
    else {Keep=false; return;}
    string SVType=getSVType(SignatureCluster);
    tuple<int,int> Site=analyzeSignatureCluster(SignatureCluster);
    int SVLen=get<1>(Site);
    Pos=get<0>(Site);//0-bsed now, after ref and alt then transform to 1-based, but should be the base before variantion. End should be the last base, but also should be transform to 1-based. So they don't change.
    //VCF version 4.2 says if alt is <ID>, the pos is the base preceding the polymorphism. No mention of the "base 1" rule.
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
}