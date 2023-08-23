import sys
import datetime

RefIndex=open(sys.argv[1],"r")
Reference=sys.argv[1][:-4]
Contigs=""
for l in RefIndex:
    sl=l.strip().split()
    Contigs+="##contig=<ID=%s,length=%s>\n"%(sl[0],sl[1])

Header="""##fileformat=VCFv4.2
##fileDate=%s
##reference=%s
%s##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">
##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number for duplications">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	bulk"""%(datetime.datetime.now().strftime("%Y%m%d"),Reference,Contigs)

SVs=[]
for BedFileName in sys.argv[3:] :
    BedFile=open(BedFileName,"r")

    for line in BedFile:
        sl=line.strip().split("\t")
        SVs.append([sl[0],int(sl[1]),int(sl[2]),sl[3],sl[4],sl[6],sl[7]])
        if SVs[-1][3]=="insertion":
            SVs[-1][1]+=1#corresponding to gensvbed.py
SVs.sort()

AllSVTypes=["DEL","DUP","INS","INV"]
if sys.argv[2]!="ALL":
    AllSVTypes=[sys.argv[2]]

print(Header)
for SV in SVs:
    ALT=None
    SVType=None
    Other=""
    Pos=SV[1]
    if (SV[3]=="insertion"):
        SVType="INS"
        ALT=SV[4]
    elif SV[3]=="deletion":
        SVType="DEL"
    elif SV[3]=="inversion":
        SVType="INV"
    elif "duplication" in SV[3]:
        SVType="DUP"
        Other=";CN="+SV[4]
    if ALT==None:
        ALT="<%s>"%SVType
    if SVType not in AllSVTypes:
        continue
    if SVType=="INS":
        SVLEN=len(ALT)
    else:
        SVLEN=int(SV[2])-int(SV[1])+1
    if SVType=="DEL":
        SVLEN=-SVLEN
    print("%s\t%s\t%s\t.\t%s\t.\tPASS\tEND=%s;SVTYPE=%s;SVLEN=%s;PRECISE%s\tGT\t%s"%(SV[0],Pos,SV[6],ALT,SV[2],SVType,SVLEN,Other,SV[5]))

BedFile.close()