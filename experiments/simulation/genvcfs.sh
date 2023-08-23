python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai ALL AllSVs.bed | bgzip -c > AllSVs.vcf.gz
tabix -f AllSVs.vcf.gz
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai DEL AllSVs.bed | bgzip -c > DELs.vcf.gz
tabix -f DELs.vcf.gz
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai INS AllSVs.bed | bgzip -c > INSs.vcf.gz
tabix -f INSs.vcf.gz
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai DUP AllSVs.bed | bgzip -c > DUPs.vcf.gz
tabix -f DUPs.vcf.gz
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai INV AllSVs.bed | bgzip -c > INVs.vcf.gz
tabix -f INVs.vcf.gz
python3 unphase.py AllSVs.bed > AllSVs.unphased.bed
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai ALL AllSVs.unphased.bed | bgzip -c > AllSVs.unphased.vcf.gz
tabix -f AllSVs.unphased.vcf.gz
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai DEL AllSVs.unphased.bed | bgzip -c > DELs.unphased.vcf.gz
tabix -f DELs.unphased.vcf.gz
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai INS AllSVs.unphased.bed | bgzip -c > INSs.unphased.vcf.gz
tabix -f INSs.unphased.vcf.gz
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai DUP AllSVs.unphased.bed | bgzip -c > DUPs.unphased.vcf.gz
tabix -f DUPs.unphased.vcf.gz
python3 genVcfFromBed.py GRCh38_chromosomes.fa.fai INV AllSVs.unphased.bed | bgzip -c > INVs.unphased.vcf.gz
tabix -f INVs.unphased.vcf.gz