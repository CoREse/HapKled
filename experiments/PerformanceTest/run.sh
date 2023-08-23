mkdir data

Ts="16 8 4 2 1"
for T in $Ts
do
    ##########################
    #ONT
    ##########################

    Reference=$DATAFOLDER/hs37d5.fa.gz
    Programs="cuteSV sniffles kled"
    if [ "$T" == "1" ];then
        Programs=$Programs" svim"
    Data=ONT
    Covs="30x"
    Ss=(5)
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
            bash runcommand.sh $ProgramName $NAME $Reference $Alignment $T $s $Others
        done
        Si=$(($Si+1))
    done

    ##########################
    #CLR
    ##########################

    Reference=$DATAFOLDER/hs37d5.fa.gz
    Programs="cuteSV sniffles kled"
    if [ "$T" == "1" ];then
        Programs=$Programs" svim"
    Data=CLR
    Covs="30x"
    Ss=(4)
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
            bash runcommand.sh $ProgramName $NAME $Reference $Alignment $T $s $Others
        done
        Si=$(($Si+1))
    done

    ##########################
    #CCS
    ##########################

    Reference=$DATAFOLDER/hs37d5.fa.gz
    Programs="cuteSV sniffles kled"
    if [ "$T" == "1" ];then
        Programs=$Programs" svim"
    Data=CCS
    Covs="30x"
    Ss=(3)
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
            bash runcommand.sh $ProgramName $NAME $Reference $Alignment $T $s $Others
        done
        Si=$(($Si+1))
    done

done