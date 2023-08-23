for var in "$@"
do
	bash getSampleSNP.sh HG00403 $var
	#echo $var
done
bcftools concat -O b HG00403*.bcf > HG00403.GRCh38.ALL.SNPs.bcf
bcftools index HG00403.GRCh38.ALL.SNPs.bcf
bcftools sort -O b HG00403.GRCh38.ALL.SNPs.bcf > HG00403.GRCh38.ALL.SNPs.sorted.bcf
bcftools index HG00403.GRCh38.ALL.SNPs.sorted.bcf
