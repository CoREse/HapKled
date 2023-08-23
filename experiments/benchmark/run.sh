Threads=8
mkdir data

##########################
#ONT
##########################

Reference=$DATAFOLDER/hs37d5.fa.gz
Programs="cuteSV sniffles svim kled"
Data=ONT
Covs="5x 10x 20x 30x 52x"
Ss=(2 3 4 5 10)
Si=0

for Cov in $Covs
do
    s=${Ss[Si]}
    NAME="HG2_"${Data}"_"${Cov}
    if [ "$Cov" == "52x" ];then
        Cov=""
    else
        Cov=.$Cov
    fi
    Alignment=$BAMFOLDER/HG002_GRCh37_ONT-UL_UCSC_20200508.phased${Cov}.bam
    for ProgramName in $Programs
    do
        bash runcommand.sh $ProgramName $NAME $Reference $Alignment $Threads $s $Others

        GREPSTR='INV\|BND\|DUP'
        SUFFIX=
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|INS'
        SUFFIX=".DEL"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.DEL.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.DEL.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|DEL'
        SUFFIX=".INS"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.INS.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.INS.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
    done
    Si=$(($Si+1))
done

##########################
#CLR
##########################

Reference=$DATAFOLDER/hs37d5.fa.gz
Programs="cuteSV sniffles svim kled"
Data=CLR
Covs="5x 10x 20x 30x 42x"
Ss=(2 3 4 4 5)
Si=0

for Cov in $Covs
do
    s=${Ss[Si]}
    NAME="HG2_"${Data}"_"${Cov}
    if [ "$Cov" == "42x" ];then
        Cov=""
    else
        Cov=.$Cov
    fi
    Alignment=$BAMFOLDER/HG002_CLR${Cov}.bam
    for ProgramName in $Programs
    do
        bash runcommand.sh $ProgramName $NAME $Reference $Alignment $Threads $s $Others

        GREPSTR='INV\|BND\|DUP'
        SUFFIX=
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|INS'
        SUFFIX=".DEL"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.DEL.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.DEL.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|DEL'
        SUFFIX=".INS"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.INS.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.INS.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
    done
    Si=$(($Si+1))
done

##########################
#CCS
##########################

Reference=$DATAFOLDER/hs37d5.fa.gz
Programs="cuteSV sniffles svim kled"
Data=CCS
Covs="5x 10x 20x 30x"
Ss=(1 2 2 3)
Si=0

for Cov in $Covs
do
    s=${Ss[Si]}
    NAME="HG2_"${Data}"_"${Cov}
    if [ "$Cov" == "30x" ];then
        Cov=""
    else
        Cov=.$Cov
    fi
    Alignment=$BAMFOLDER/HG002_CCS$Cov.bam
    for ProgramName in $Programs
    do
        bash runcommand.sh $ProgramName $NAME $Reference $Alignment $Threads $s $Others

        GREPSTR='INV\|BND\|DUP'
        SUFFIX=
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|INS'
        SUFFIX=".DEL"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.DEL.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.DEL.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|DEL'
        SUFFIX=".INS"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.INS.vcf.gz"
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.INS.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
    done
    Si=$(($Si+1))
done
##########################
#SIMs
##########################

Reference=../simulation/GRCh38_chromosomes.fa
Programs="cuteSV sniffles svim kled"
Data=SIM
Covs="5x 10x 20x 30x"
Ss=(2 3 4 4)
Si=0

Chr="chr1"

for Cov in $Covs
do
    s=${Ss[Si]}
    NAME="HG2_"${Data}${Chr}"_"${Cov}
    Cov=_$Cov
    if [ "$Chr" == "chr1" ];then
        Reference=../simulation/GRCh38_chr1.fa
    fi
    Alignment=../simulation/SG${Chr}_bam${Cov}/sim.srt${Detag}.bam
    for ProgramName in $Programs
    do
        bash runcommand.sh $ProgramName $NAME $Reference $Alignment $Threads $s $OtherParas
        
        if [ "$Chr" == "chr1" ];then
            INCLUDEBED=" --includebed ../simulation/GRCh38_chromosomes.1.bed"
        fi
        GREPSTR='BND'
        SUFFIX=
        GS=../simulation/AllSVs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|INS'
        SUFFIX=".DEL"
        GS=../simulation/DELs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|DEL'
        SUFFIX=".INS"
        GS=../simulation/INSs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INS\|INV\|BND\|DEL'
        SUFFIX=".DUP"
        GS=../simulation/DUPs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INS\|DUP\|BND\|DEL'
        SUFFIX=".INV"
        GS=../simulation/INVs.unphased.vcf.gz$INCLUDEBED
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
    done
    Si=$(($Si+1))
done