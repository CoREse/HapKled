#!/usr/bin/bash
Threads=32
mkdir data

if [ -z "$DataDir" ];then
    >&2 echo "`date '+[%Y-%m-%d %H:%M:%S]'` You should first provide the path of the data directory that contains all the required data by setting the environment variable DataDir."
fi
##########################
#ONT
##########################
GSG=""

DataSet="ONT"
Programs="HapKled kled cuteSV sniffles duet"
Covs="30x 20x 10x 5x"
Ss=(5 4 3 2)
Si=0

for Cov in $Covs
do
    s=${Ss[Si]}
    NAME="HG2_"${DataSet}${Chr}${Detag}"_"${Cov}
    DataType="ONT"
    Reference=$DataDir/hs37d5.fa
    Alignment="$DataDir/bams/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.${Cov}.bam"
    INCLUDEBED=" --includebed $DataDir/GIAB/HG002_SVs_Tier1_v0.6.bed"
    Cov=_$Cov
    for ProgramName in $Programs
    do
        bash runcommand.sh $ProgramName $DataType $NAME $Reference $Alignment $Threads $s $OtherParas
    done
    for ProgramName in $Programs
    do
        GREPSTR='=BND\|=INV\|=DUP'
        SUFFIX=
        GS=$DataDir/GIAB/HG002_SVs_Tier1_v0.6.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $Reference $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='=INV\|=BND\|=DUP\|=INS'
        SUFFIX=".DEL"
        GS=$DataDir/GIAB/HG002_SVs_Tier1_v0.6.DEL.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $Reference $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='=INV\|=BND\|=DUP\|=DEL'
        SUFFIX=".INS"
        GS=$DataDir/GIAB/HG002_SVs_Tier1_v0.6.INS.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $Reference $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
    done
    Si=$(($Si+1))
done

##########################
#SIMONT
##########################


DataSet="SIMONT"
Programs="HapKled kled cuteSV sniffles duet"
Covs="30x"
Ss=(5)

Si=0

for Cov in $Covs
do
    s=${Ss[Si]}
    NAME="HG2_"${DataSet}${Chr}${Detag}"_"${Cov}
    DataType="ONT"
    Reference=$DataDir/GRCh38_chromosomes.fa
    if [ "$Cov" != "" ];then
        Alignment="$DataDir/bams/SG_bam_${Cov}/sim.srt.bam"
    else
        continue
    fi
    INCLUDEBED=""
    Cov=_$Cov
    for ProgramName in $Programs
    do
        bash runcommand.sh $ProgramName $DataType $NAME $Reference $Alignment $Threads $s $OtherParas
    done
    for ProgramName in $Programs
    do
        GREPSTR='=BND'
        SUFFIX=
        GS=/home/user/zhangzhendong/workspace/kled/exp/wholehybrid/starredsucceededrunPbWhole4171692cu26/work/AllSVs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $Reference $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='=INV\|=BND\|=DUP\|=INS'
        SUFFIX=".DEL"
        GS=/home/user/zhangzhendong/workspace/kled/exp/wholehybrid/starredsucceededrunPbWhole4171692cu26/work/DELs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $Reference $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='=INV\|=BND\|=DUP\|=DEL'
        SUFFIX=".INS"
        GS=/home/user/zhangzhendong/workspace/kled/exp/wholehybrid/starredsucceededrunPbWhole4171692cu26/work/INSs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $Reference $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='=INS\|=INV\|=BND\|=DEL'
        SUFFIX=".DUP"
        GS=/home/user/zhangzhendong/workspace/kled/exp/wholehybrid/starredsucceededrunPbWhole4171692cu26/work/DUPs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $Reference $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='=INS\|=DUP\|=BND\|=DEL'
        SUFFIX=".INV"
        GS=/home/user/zhangzhendong/workspace/kled/exp/wholehybrid/starredsucceededrunPbWhole4171692cu26/work/INVs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $Reference $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
    done
    Si=$(($Si+1))
done