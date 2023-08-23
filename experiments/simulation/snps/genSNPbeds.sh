bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' $1 | grep "1|0\|1|1" | awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' > HACk.snps.h1.bed
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' $1 | grep "0|1\|1|1" | awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' > HACk.snps.h2.bed
