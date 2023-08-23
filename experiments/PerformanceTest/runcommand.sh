ProgramName=$1
NAME=$2
Reference=$3
Alignment=$4
Threads=$5
s=$6
Others=$7

sudo clearcache
if [ "$ProgramName" = "kled" ]
then
  KledParas=
  if [[ $Alignment =~ "ONT" ]]; then
    KledParas=""
  elif [[ $Alignment =~ "CLR" ]]; then
    KledParas="--CLR"
  elif [[ $Alignment =~ "CCS" ]]; then
    KledParas="--CCS"
  else
    KledParas=""
  fi
tmemc -d -p ${ProgramName} &
(set -x;\time -v ${ProgramName} -R $Reference $Alignment -t $Threads $KledParas > data/${ProgramName}_$NAME.vcf)
elif [ "$ProgramName" = "cuteSV" ]
then
  mkdir cuteWork
  rm -r cuteWork/*
  CuteParas=""
  # For PacBio CLR data:
	# 	--max_cluster_bias_INS		100
	# 	--diff_ratio_merging_INS	0.3
	# 	--max_cluster_bias_DEL	200
	# 	--diff_ratio_merging_DEL	0.5

	# For PacBio CCS(HIFI) data:
	# 	--max_cluster_bias_INS		1000
	# 	--diff_ratio_merging_INS	0.9
	# 	--max_cluster_bias_DEL	1000
	# 	--diff_ratio_merging_DEL	0.5

	# For ONT data:
	# 	--max_cluster_bias_INS		100
	# 	--diff_ratio_merging_INS	0.3
	# 	--max_cluster_bias_DEL	100
	# 	--diff_ratio_merging_DEL	0.3
  if [[ $Reference == *.gz ]]; then
    # Remove ".gz" from filename since cuteSV doesn't support .gz reference files
    Reference=${Reference::-3}
  fi
  if [[ $Alignment =~ "ONT" ]]; then
    CuteParas=""
  elif [[ $Alignment =~ "CLR" ]]; then
    CuteParas="--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 -md 500 -mi 500"
  elif [[ $Alignment =~ "CCS" ]]; then
    CuteParas="--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -md 500 -mi 500"
  else
    CuteParas="--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -md 0 -mi 0"
  fi
  tmemc -d -p $ProgramName &
  (set -x;\time -v $ProgramName $Alignment $Reference data/${ProgramName}_$NAME.vcf --genotype ./cuteWork/ -s $s $CuteParas -t $Threads $Others)
elif [ "$ProgramName" = "sniffles" ]
then
tmemc -d -p $ProgramName &
(set -x;\time -v ${ProgramName} --input $Alignment --vcf data/${ProgramName}_$NAME.vcf --threads $Threads --allow-overwrite)
elif [ "$ProgramName" = "NanoSV" ]
then
tmemc -d -p $ProgramName &
(set -x;\time -v NanoSV $Alignment -c nanoconfig.txt -t $Threads > data/${ProgramName}_$NAME.vcf)
rm data/${ProgramName}_$NAME.unsorted.vcf
elif [ "$ProgramName" = "svim" ]
then
tmemc -d -p $ProgramName &
mkdir svimWork
rm -r svimWork/*
(set -x;\time -v svim alignment svimWork $Alignment $Reference > data/${ProgramName}_$NAME.vcf)
cp svimWork/variants.vcf data/${ProgramName}_$NAME.vcf
fi