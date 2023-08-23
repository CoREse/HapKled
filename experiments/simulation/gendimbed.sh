Ref=$1
cut -f1,2 $2/*.fai $Ref.fai > $2.haplochroms.dim.tsv
cat $2.haplochroms.dim.tsv | sort  | awk '$2 > maxvals[$1] {lines[$1]=$0; maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' > $2.maxdims.tsv
awk 'OFS=FS="\t"''{print $1, "1", $2, "100.0", "100.0"}' $2.maxdims.tsv > $2.shorts.laser.simple.bed