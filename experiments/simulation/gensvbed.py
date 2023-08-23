import sys
import random
import pysam

ShortInsMean=45
ShortInsStd=22
ShortInsRatio=0.5
MediumInsMean=350
MediuInsStd=200
MediuInsRatio=0.4
LargeInsMin=1000
LargeInsMax=100000

CHRPrefix="chr"

Alphabet="ATGC"
def genInsBases(SVLen):
    Base=""
    for i in range(SVLen):
        Base+=Alphabet[random.randrange(4)]
    return Base

def genAF(Min=0.01):
    AF=random.gammavariate(0.2,1)
    if AF<Min:
        AF=Min
    if AF>1:
        AF=1
    return AF

class SV:
    def __init__(self,Type="", Chr="",Start=0, End=0, SVLen=0, AF=0, ID=0):
        self.Type=Type
        self.Start=Start
        self.End=End
        self.Chr=Chr
        self.AF=AF
        self.ID=ID
        if SVLen==0:
            self.SVLen=End-Start
        else:
            self.SVLen=SVLen
        if (Type=="insertion"):
            self.genIns()
        elif (Type=="tandem duplication"):
            self.genCN()
    def genIns(self):
        if self.End==self.Start and self.SVLen==0:
            Rnd=random.uniform(0,1)
            if (Rnd<=ShortInsRatio):
                self.SVLen=random.gauss(ShortInsMean,ShortInsStd)
            elif (Rnd<=MediuInsRatio):
                self.SVLen=random.gauss(MediumInsMean,MediuInsStd)
            else:
                self.SVLen=random.uniform(LargeInsMin,LargeInsMax)
            self.SVLen=int(self.SVLen)
        self.InsBases=genInsBases(self.SVLen)
    def genCN(self):
        CN=random.gauss(2,2)
        CN=int(CN)
        if (CN<2):
            CN=2
        self.CN=CN
        #The SVLen is from insertion, so to devide by CN
        self.SVLen/=CN-1
        if self.SVLen<30:
            self.SVLen+=30
        self.SVLen=int(self.SVLen)
        self.End=self.Start+self.SVLen
    def __lt__(self,other):
        if (self.Chr==other.Chr):
            return self.Start<other.Start
        try:
            ChrN=int(self.Chr)
            OChrN=int(other.Chr)
            return ChrN<OChrN
        except:
            return self.Chr<other.Chr
    def __eq__(self,other):
        return self.ID==other.ID
    def __str__(self):
        Str="%s\t"%(self.Chr if len(self.Chr)>2 else CHRPrefix+self.Chr)
        if (self.Type=="insertion"):
            Str+="%s\t%s\t"%(self.Start-1,self.Start)
        else:
            Str+="%s\t%s\t"%(self.Start,self.End)
        Str+="%s\t"%(self.Type)
        if (self.Type=="tandem duplication"):
            Str+="%s\t"%(self.CN)
        elif (self.Type=="insertion"):
            Str+=self.InsBases+"\t"
        else:
            Str+="None\t"
        Str+="0"
        return Str

def getSVLenDist(VCFFileName):
    SVLenDist={"DEL":[],"INS":[]}
    vcf=pysam.VariantFile(VCFFileName,"r")
    for record in vcf.fetch():
        SVTYPE=record.info["SVTYPE"]
        for SVLEN in record.info["SVLEN"]:
            if (SVTYPE=="DEL"):
                SVLenDist["DEL"].append(abs(int(SVLEN)))
            elif (SVTYPE=="INS"):
                SVLenDist["INS"].append(abs(int(SVLEN)))
            else:
                continue
    vcf.close()
    SVLenDist["DEL"]=sorted(SVLenDist["DEL"])
    SVLenDist["INS"]=sorted(SVLenDist["INS"])
    return SVLenDist

def deOverlap(SVs):
    SVs.sort()
    i=0
    while i< len(SVs)-1:
        if (i+1<len(SVs) and SVs[i].Chr==SVs[i+1].Chr and SVs[i].End>=SVs[i+1].Start):
            del SVs[i+1]
            continue
        i+=1

def genHaplotypes(SVs,H1SVs,H2SVs):
    for V in SVs:
        P=V.AF*V.AF/(1.0-(1.0-V.AF)**2)
        if random.uniform(0,1)<P:#homo
            H1SVs.append(V)
            H2SVs.append(V)
        else:
            if random.uniform(0,1)<0.5:
                H1SVs.append(V)
            else:
                H2SVs.append(V)

def removeDuplicants(SVs):
    SVs=SVs.copy()
    SVs.sort()
    i=0
    while i< len(SVs)-1:
        if (i+1<len(SVs) and SVs[i]==SVs[i+1]):
            del SVs[i+1]
            continue
        i+=1  
    return SVs

def getChrs(FaFileName):
    Chroms=[]
    TotalLength=0
    fa=pysam.FastaFile(FaFileName)
    for chrom in fa.references:
        Chroms.append((chrom,fa.get_reference_length(chrom)))
        TotalLength+=fa.get_reference_length(chrom)
    return Chroms,TotalLength

def getPosChr(Pos,Chroms):
    for i in range(len(Chroms)):
        if Pos<Chroms[i][1]:
            return Chroms[i][0],Pos
        Pos-=Chroms[i][1]
    return -1,-1

def find(SortedArray, Value):
    for i in range(len(SortedArray)):
        if Value<=SortedArray[i]:
            return i
    return len(SortedArray)

TypeNames=["DEL","INS","DUP","INV"]
TypeLongNames=["deletion","insertion","tandem duplication","inversion"]
def genSVs(Chroms,TotalLength, SVLenDist,Numbers=(70000,70000,5000,5000)):
    SVs=[]
    for Typei in range(4):
        Type=TypeNames[Typei]
        for i in range(Numbers[Typei]):
            RandI=random.randint(0,len(SVLenDist[Type])-1)
            SVLEN=SVLenDist[Type][RandI]
            if RandI!=len(SVLenDist[Type])-1:
                RandR=random.random()
                SVLEN=SVLEN*RandR+SVLenDist[Type][RandI+1]*(1-RandR)
            SVLEN=SVLEN+((-1)**random.randint(0,1))*SVLEN*random.random()*0.1
            SVLEN=int(SVLEN)
            if SVLEN<30:
                SVLEN=30
            AF=genAF()
            GPos=random.randint(0,TotalLength-1)
            Chr,Pos=getPosChr(GPos,Chroms)
            if Chr==-1:
                raise -1
            # print(Chr,Pos)
            SVs.append(SV(TypeLongNames[Typei],Chr,Pos,Pos+SVLEN,SVLEN,AF,"SIMULATED.%s.%s"%(Type,i)))
    return SVs

import sys
def run():
    TotalSVN=50000
    Chroms, TotalLength=getChrs(sys.argv[1])
    SVLenDist=getSVLenDist(sys.argv[2])
    H1BedName=sys.argv[3]
    H2BedName=sys.argv[4]
    Seed=int(sys.argv[5])
    random.seed(Seed)
    print("Seed: %s."%Seed, file=sys.stderr)
    MinDupLen=500
    MinInvLen=200
    SVLenDist["DUP"]=SVLenDist["INS"][find(SVLenDist["INS"],MinDupLen):]
    SVLenDist["INV"]=SVLenDist["DEL"][find(SVLenDist["DEL"],MinInvLen):]
    SVs=genSVs(Chroms,TotalLength,SVLenDist)
    SVs=random.sample(SVs,TotalSVN)
    H1SVs=[]
    H2SVs=[]
    genHaplotypes(SVs,H1SVs,H2SVs)
    deOverlap(H1SVs)
    deOverlap(H2SVs)
    AllSVs=removeDuplicants(H1SVs+H2SVs)

    NDEL=0
    NINS=0
    NDUP=0
    NINV=0

    for V in AllSVs:
        if V.Type=="deletion":
            NDEL+=1
        elif V.Type=="insertion":
            NINS+=1
        elif V.Type=="tandem duplication":
            NDUP+=1
        elif V.Type=="inversion":
            NINV+=1

    print("All SVs: %d, SVs for H1: %d, SVs for H2: %d. DEL, INS, DUP and INV: %d %d %d %d"%(len(AllSVs),len(H1SVs),len(H2SVs),NDEL,NINS,NDUP,NINV),file=sys.stderr)

    H1SVIDs=set()
    H2SVIDs=set()

    H1F=open(H1BedName,"w")
    for V in H1SVs:
        print(V,file=H1F)
        H1SVIDs.add(V.ID)
    H1F.close()
    H2F=open(H2BedName,"w")
    for V in H2SVs:
        print(V,file=H2F)
        H2SVIDs.add(V.ID)
    H2F.close()

    for V in AllSVs:
        GT=""
        if (V.ID in H1SVIDs):
            GT+="1"
        else:
            GT+="0"
        GT+="|"
        if (V.ID in H2SVIDs):
            GT+="1"
        else:
            GT+="0"
        print(str(V)+"\t"+GT+"\t"+V.ID)

if __name__=="__main__":
    run()
