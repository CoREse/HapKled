IN=$1
OUT=$2
DISCORDANT_VARIANTS=$3
jasmine --output_genotypes --pre_normalize file_list=$IN out_file=$OUT
python3 calcMDR.py $OUT $DISCORDANT_VARIANTS