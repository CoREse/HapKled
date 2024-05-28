ProgramName=$1
Ref=$2
GREPSTR=$3
GS=$4
NAME=$5
s=$6
SUFFIX=$7

PASSONLY=" --passonly"

grep -v ${GREPSTR} $NAME.vcf | bgzip -c > $NAME$SUFFIX.vcf.gz

tabix -f $NAME$SUFFIX.vcf.gz
rm -r $NAME$SUFFIX.cmp

VERI=truvari
(set -x;$VERI bench -b $GS -c $NAME$SUFFIX.vcf.gz --reference $Ref -p 0.00$PASSONLY -o $NAME$SUFFIX.cmp)
