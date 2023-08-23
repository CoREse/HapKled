mv HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam HG002_CCS.bam
samtools index -@7 HG002_CCS.bam
mv HG002_PacBio_GRCh37.bam HG002_CLR.bam
samtools index -@7 HG002_CLR.bam
samtools view -s 0.33 -b -h -@ 7 HG002_CCS.bam > HG002_CCS.10x.bam && samtools index -@ 7 HG002_CCS.10x.bam &
samtools view -s 0.66 -b -h -@ 7 HG002_CCS.bam > HG002_CCS.20x.bam && samtools index -@ 7 HG002_CCS.20x.bam &
samtools view -s 0.165 -b -h -@ 7 HG002_CCS.bam > HG002_CCS.5x.bam && samtools index -@ 7 HG002_CCS.5x.bam &
samtools view -s 0.165 -b -h -@ 7 HG002_CLR.bam > HG002_CLR.10x.bam && samtools index -@ 7 HG002_CLR.10x.bam &
samtools view -s 0.33 -b -h -@ 7 HG002_CLR.bam > HG002_CLR.20x.bam && samtools index -@7 HG002_CLR.20x.bam &
samtools view -s 0.5 -b -h -@ 7 HG002_CLR.bam > HG002_CLR.30x.bam && samtools index -@ 7 HG002_CLR.30x.bam &
samtools view -s 0.0825 -b -h -@ 7 HG002_CLR.bam > HG002_CLR.5x.bam &
samtools view -s 0.1923 -b -h -@ 7 HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam > HG002_GRCh37_ONT-UL_UCSC_20200508.phased.10x.bam
samtools view -s 0.3846 -b -h -@ 7 HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam > HG002_GRCh37_ONT-UL_UCSC_20200508.phased.20x.bam
samtools view -s 0.576923 -b -h -@ 7 HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam > HG002_GRCh37_ONT-UL_UCSC_20200508.phased.30x.bam
samtools view -s 0.09615 -b -h -@ 7 HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam > HG002_GRCh37_ONT-UL_UCSC_20200508.phased.5x.bam