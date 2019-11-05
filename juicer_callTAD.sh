#!/bin/bash
juicertool="java -Xms512m -Xmx2048m -jar $(cd $(dirname $0) && pwd)/../binaries/Aidenlab/juicer_tools.1.8.9_jcuda.0.8.jar"

norm=$1
hic=$2
odir=$3

dir=$odir/TAD/$norm
mkdir -p $dir

for res in 10000 25000 50000; do
    $juicertool arrowhead --ignore_sparsity -m 2000 -r $res -k $norm $hic $dir
    #	cat $dir/${res}_blocks.bedpe | awk '{OFS="\t"} NR>1 {print $1, $2, $3 }' | sort -k1 -k2,2n > tmp1.txt
    #	perl /work/sakata/Perl/border.pl tmp1.txt $res | sort -k1 -k2,2n | uniq > tmp2.txt
    #	addchr.pl tmp2.txt > $dir/borders_${res}_blocks
    #	addchr.pl $dir/${res}_blocks.bedpe | awk '{OFS="\t"} NR>1 {print }' > $dir/${res}_blocks.txt
#	rm -f tmp*.txt
done
