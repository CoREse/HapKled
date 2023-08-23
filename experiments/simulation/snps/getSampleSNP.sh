Name=$2
chr=(${Name//./ })
chr=${chr[1]}
#echo $chr
echo Getting $chr for $1
(set -x;bcftools view -O b -o $1.$chr.SNPs.bcf -s $1 -m2 -M2 -c 1 -C 2 --threads 8 $2)
bcftools index $1.$chr.SNPs.bcf
