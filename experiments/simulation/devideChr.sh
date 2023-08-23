FullGenome=$1
Chr=$2
mkdir $FullGenome$Chr
samtools faidx $FullGenome/h1.fa $Chr > $FullGenome$Chr/h1.fa
samtools faidx $FullGenome$Chr/h1.fa
samtools faidx $FullGenome/h2.fa $Chr > $FullGenome$Chr/h2.fa
samtools faidx $FullGenome$Chr/h2.fa