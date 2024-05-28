# Kled: a haplotype-aware structural variant calling approach for Oxford Nanopore sequencing data
## Usage
Before running HapKled you should first compile (see compiling the haplotype-aware kled) the haplotype-aware kled and put the path of it to the environment variable HapAwareKled.
```
export HapAwareKled=/path/to/hap-aware-kled
```
And you also need an installed Clair3 and Whatshap, and export the path of the Clair3 models to environment variable Clair3ModelPath.
```
export Clair3ModelPath=/path/to/bin/models
```
HapKled need a reference file (fasta) and at least one bam (sam/bam/cram) file that stores the mapped reads to call SVs, and output a VCF file to the standard output.
```
HapKled -R Refernce.fa Sample.bam > SVs.vcf
```

For the description of all parameters:
```
HapKled --help
```

## Compiling the haplotype-aware kled
Dependencies: openmp and dependencies of htslib (-lz -lm -lbz2 -llzma -lcurl -lpthread -lcrypto -ldeflate)

Here are instructions to get some dependencies from source if systemwide installation is not available:
- zlib: download from https://github.com/madler/zlib/releases, unzip and cd and ./configure && make && make install prefix=/path/to/your/local
- libbz2: download from https://sourceforge.net/projects/bzip2/, unzip and cd and make && make install PREFIX=/path/to/your/local
And export CPATH to make compiler know where the lib is: export CPATH=/path/to/your/local/include:$CPATH

You also need a compiler that supports C++ 17 standard.

To get the kled binary:
```
#to get htslib
git submodule update --init --recursive

#make kled
make
```

If you want to install kled:
```
make install
#or if you want to install to a place other than /usr/local:
PREFIX=PATH_YOU_SELECT make install
```

## Experiments
You should run the run.sh in the experiments/simulation/run.sh to get the simulated bam, and along with it, put the downloaded bams and fastas to the $DataDir.

Then run experiments/benchmark/run.sh to get the benchmark results.