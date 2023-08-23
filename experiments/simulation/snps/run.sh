wget -i snpurls.txt
bash getAll.sh *.vcf.gz
bash genSNPbeds.sh HG00403.GRCh38.ALL.SNPs.sorted.bcf
