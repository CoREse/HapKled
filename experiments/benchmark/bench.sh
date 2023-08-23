ProgramName=$1
GREPSTR=$2
GS=$3
NAME=$4
s=$5
SUFFIX=$6

sudo clearcache

if [ "$ProgramName" = "cuteSV" ]
then
    GREPSTR=cuteSV.$GREPSTR
    GREPSTR=`echo $GREPSTR | sed -e "s/|/|cuteSV./g"`
    GREPSTR=$GREPSTR'\|STRAND=--'
fi
PASSONLY=" --passonly"

if [ "$ProgramName" = "NanoSV" ] || [ "$ProgramName" = "svim" ]
then
if [ ! -f $NAME.unsorted.vcf ];then
    mv $NAME.vcf $NAME.unsorted.vcf
    bcftools sort $NAME.unsorted.vcf > $NAME.vcf
else
    mv $NAME.vcf $NAME.temp.vcf
    bcftools sort $NAME.temp.vcf > $NAME.vcf
fi
fi

if [ "$ProgramName" = "svim" ]
then
    SUP=$s
    echo "svim sup: ${SUP}"
    bcftools filter -i "INFO/SUPPORT>=${SUP}" $NAME.vcf | grep -v ${GREPSTR} | python3 addsvlenforsviminv.py | python3 altdupforsvim.py | bgzip -c > $NAME$SUFFIX.vcf.gz
elif [ "$ProgramName" = "NanoSV" ];then
    grep -v ${GREPSTR} $NAME.vcf | sed "s/##INFO=<ID=RT,Number=3/##INFO=<ID=RT,Number=1/g" | bgzip -c > $NAME$SUFFIX.vcf.gz
else
    grep -v ${GREPSTR} $NAME.vcf | bgzip -c > $NAME$SUFFIX.vcf.gz
fi

tabix -f $NAME$SUFFIX.vcf.gz
rm -r $NAME$SUFFIX.cmp
#truvari version: Truvari v2.1 - Structural Variant Benchmarking and Annotation
Ref=$DATAFOLDER/hs37d5.fa.gz
if [[ $NAME =~ "SIM" ]]; then
    Ref=../simulation/GRCh38_chromosomes.fa
fi
VERI=truvari
(set -x;$VERI bench -b $GS -c $NAME$SUFFIX.vcf.gz --reference $Ref -r 1000 -p 0.00$PASSONLY -o $NAME$SUFFIX.cmp)
