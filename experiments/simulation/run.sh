Seed=1000
SimulatedName="SG"
bash prepareFastas.sh
Ref=GRCh38_chromosomes.fa
cd snps
bash run.sh 2>run.log 1>&2
cd ..
bash genSVBeds.sh $Seed
bash genvcfs.sh
bash genfa.sh $SimulatedName
bash devideChr.sh $SimulatedName chr1
Ref=GRCh38_chr1.fa
SimulatedName="SGchr1"
bash gendimbed.sh $Ref $SimulatedName
bash genreads.sh $Ref $SimulatedName