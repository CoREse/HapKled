ProgramName=$1
DataType=$2
NAME=$3
Reference=$4
Alignment=$5
Threads=$6
s=$7
Others=$8

if [ "$ProgramName" = "kled" ]
then
  KledParas=
  if [[ $DataType =~ "ONT" ]]; then
    KledParas=""
  elif [[ $DataType =~ "CLR" ]]; then
    KledParas="--CLR"
  elif [[ $DataType =~ "CCS" ]]; then
    KledParas="--CCS"
  else
    KledParas=""
  fi
  (set -x;\time -v kled -R $Reference $Alignment -t $Threads $KledParas > data/${ProgramName}_$NAME.vcf)
elif [ "$ProgramName" = "cuteSV" ]
then
  mkdir cuteWork
  rm -r cuteWork/*
  if [[ $Reference == *.gz ]]; then
    Reference=${Reference::-3}
  fi
  CuteParas=""
  if [[ $DataType =~ "ONT" ]]; then
    CuteParas=""
  elif [[ $DataType =~ "CLR" ]]; then
    CuteParas="--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 -md 500 -mi 500"
  elif [[ $DataType =~ "CCS" ]]; then
    CuteParas="--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -md 500 -mi 500"
  fi
  (set -x;\time -v cuteSV $Alignment $Reference data/${ProgramName}_$NAME.vcf --genotype ./cuteWork/ -s $s $CuteParas -t $Threads $Others)
elif [ "$ProgramName" = "sniffles" ]
then
(set -x;\time -v ${ProgramName} --input $Alignment --vcf data/${ProgramName}_$NAME.vcf --threads $Threads --allow-overwrite)
elif [ "$ProgramName" = "PhaseSV" ]
then
  KledParas=
  if [[ $DataType =~ "ONT" ]]; then
    KledParas=""
  elif [[ $DataType =~ "CLR" ]]; then
    KledParas="--callerparas --CLR"
  elif [[ $DataType =~ "CCS" ]]; then
    KledParas="--callerparas --CCS"
  else
    KledParas=""
  fi
  (set -x;\time -v HapKled -R $Reference $Alignment -t $Threads $KledParas --workdir ./HapKledWork_$NAME > data/${ProgramName}_$NAME.vcf)
elif [ "$ProgramName" = "duet" ]
then
# You should . condainit.sh here to make the conda work
# conda activate duet
mkdir duetWork
rm -r duetWork/*
duet $Alignment $Reference duetWork -t 32
sed 's/HP/GT/g' duetWork/phased_sv.vcf | sed 's/##INFO=<ID=SVLEN/##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n##INFO=<ID=END,Number=1,Type=Integer,Description="Estimated end of the variant">\n##INFO=<ID=SVLEN/' | sed 's/=<DEL>/=DEL/g' | sed 's/=<INS>/=INS/g' | sed 's/=<DUP>/=DUP/g' | sed 's/=<INV>/=INV/g' > data/${ProgramName}_$NAME.vcf
fi