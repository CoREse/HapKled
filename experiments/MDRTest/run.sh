mkdir data
##########################
#ONT
##########################

Reference=$DATAFOLDER/hs37d5.fa.gz
Programs="cuteSV sniffles svim kled"
Data=ONT
Ts="8"


for ProgramName in $Programs
do
    for T in $Ts
    do
        NAME="HG2_"${Data}"_52x"
        s=10
        bash runcommand.sh $ProgramName $NAME $Reference $BAMFOLDER/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam $T 10 $Others
        GREPSTR="BND"
        SUFFIX=
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|INS'
        SUFFIX=".DEL"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|DEL'
        SUFFIX=".INS"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INS\|INV\|BND\|DEL'
        SUFFIX=".DUP"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INS\|DUP\|BND\|DEL'
        SUFFIX=".INV"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX

        NAME="HG3_"${Data}"_90x"
        s=16
        bash runcommand.sh $ProgramName $NAME $Reference $BAMFOLDER/HG003_GRCh37_ONT-UL_UCSC_20200508.bam $T $s $Others
        GREPSTR="BND"
        SUFFIX=
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|INS'
        SUFFIX=".DEL"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|DEL'
        SUFFIX=".INS"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INS\|INV\|BND\|DEL'
        SUFFIX=".DUP"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INS\|DUP\|BND\|DEL'
        SUFFIX=".INV"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX

        NAME="HG4_"${Data}"_93x"
        s=16
        bash runcommand.sh $ProgramName $NAME $Reference $BAMFOLDER/HG004_GRCh37_ONT-UL_UCSC_20200508.bam $T $s $Others
        GREPSTR="BND"
        SUFFIX=
        GS="$DATAFOLDER/HG002_SVs_Tier1_v0.6.vcf.gz --includebed $DATAFOLDER/HG002_SVs_Tier1_v0.6.bed"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|INS'
        SUFFIX=".DEL"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INV\|BND\|DUP\|DEL'
        SUFFIX=".INS"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INS\|INV\|BND\|DEL'
        SUFFIX=".DUP"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
        GREPSTR='INS\|DUP\|BND\|DEL'
        SUFFIX=".INV"
        bash bench.sh $ProgramName $GREPSTR "$GS" data/${ProgramName}_$NAME $s $SUFFIX
    done
    SUFFIX=
    gunzip -c data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf
    echo -e "data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf" > input.txt
    (set -x;bash calcMDR.sh input.txt data/${ProgramName}_Merge${SUFFIX}.vcf data/${ProgramName}_DiscordantVariants${SUFFIX}.vcf)
    SUFFIX=.DEL
    gunzip -c data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf
    echo -e "data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf" > input.txt
    (set -x;bash calcMDR.sh input.txt data/${ProgramName}_Merge${SUFFIX}.vcf data/${ProgramName}_DiscordantVariants${SUFFIX}.vcf)
    SUFFIX=.INS
    gunzip -c data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf
    echo -e "data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf" > input.txt
    (set -x;bash calcMDR.sh input.txt data/${ProgramName}_Merge${SUFFIX}.vcf data/${ProgramName}_DiscordantVariants${SUFFIX}.vcf)
    SUFFIX=.DUP
    gunzip -c data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf
    echo -e "data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf" > input.txt
    (set -x;bash calcMDR.sh input.txt data/${ProgramName}_Merge${SUFFIX}.vcf data/${ProgramName}_DiscordantVariants${SUFFIX}.vcf)
    SUFFIX=.INV
    gunzip -c data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf
    gunzip -c data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz > data/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf
    echo -e "data/${ProgramName}_"HG2_"${Data}_52x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG3_"${Data}_90x${SUFFIX}.vcf.gz.vcf\ndata/${ProgramName}_"HG4_"${Data}_93x${SUFFIX}.vcf.gz.vcf" > input.txt
    (set -x;bash calcMDR.sh input.txt data/${ProgramName}_Merge${SUFFIX}.vcf data/${ProgramName}_DiscordantVariants${SUFFIX}.vcf)
done
