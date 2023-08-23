ProgramName=$1
GREPSTR=$2
GS=$3
NAME=$4
s=$5
SUFFIX=$6

sudo clearcache

PASSONLY=""
echo "bench $@"

if [ "$ProgramName" = "cuteSV" ]
then
GREPSTR=cuteSV.$GREPSTR
GREPSTR=`echo $GREPSTR | sed -e "s/|/|cuteSV./g"`
fi
if [ "$ProgramName" = "cuteSV" ] || [ "$ProgramName" = "sniffles" ] || [ "$ProgramName" = "svim" ]
then
PASSONLY=" --passonly"
fi

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
bcftools filter -i "INFO/SUPPORT>=${SUP}" $NAME.vcf | grep -v ${GREPSTR} | python3 ../addsvlenforsviminv.py | python3 ../altdupforsvim.py | bgzip -c > $NAME$SUFFIX.vcf.gz
elif [ "$ProgramName" = "NanoSV" ];then
grep -v ${GREPSTR} $NAME.vcf | sed "s/##INFO=<ID=RT,Number=3/##INFO=<ID=RT,Number=1/g" | bgzip -c > $NAME$SUFFIX.vcf.gz
else
grep -v ${GREPSTR} $NAME.vcf | bgzip -c > $NAME$SUFFIX.vcf.gz
fi

tabix -f $NAME$SUFFIX.vcf.gz
rm -r $NAME$SUFFIX.cmp
