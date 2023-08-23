NAME=$1
mkdir ${NAME}_bam_20x
samtools view -s 0.66666 -b -h -@ 7 ${NAME}_bam_30x/sim.srt.bam > ${NAME}_bam_20x/sim.srt.bam
samtools index ${NAME}_bam_20x/sim.srt.bam
mkdir ${NAME}_bam_10x
samtools view -s 0.33333 -b -h -@ 7 ${NAME}_bam_30x/sim.srt.bam > ${NAME}_bam_10x/sim.srt.bam
samtools index ${NAME}_bam_10x/sim.srt.bam
mkdir ${NAME}_bam_5x
samtools view -s 0.166665 -b -h -@ 7 ${NAME}_bam_30x/sim.srt.bam > ${NAME}_bam_5x/sim.srt.bam
samtools index ${NAME}_bam_5x/sim.srt.bam