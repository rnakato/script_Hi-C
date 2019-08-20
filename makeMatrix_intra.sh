#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <norm> <matrixdir> <hic file> <binsize> <build>" 1>&2
}

if [ $# -ne 5 ]; then
  usage
  exit 1
fi

norm=$1
matrixdir=$2
hic=$3
binsize=$4
build=$5

gt=/work/Database/UCSC/$build/genome_table

if test $build = "mm10" -o $build = "mm9"; then
    chrnum=19
else
    chrnum=22
fi

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx2048m -jar $pwd/../binaries/Aidenlab/juicer_tools.1.8.9_jcuda.0.8.jar"
dir=$matrixdir/intrachromosomal/$binsize
mkdir -p $dir

for chr in $(seq 1 $chrnum) X
do
    echo "chr$chr"
    for type in observed #oe
    do
        $juicertool dump $type $norm $hic $chr $chr BP $binsize $dir/$type.$norm.chr$chr.txt
        $pwd/convert_JuicerDump_to_dense.py \
	    $dir/$type.$norm.chr$chr.txt \
	    $dir/$type.$norm.chr$chr.matrix.gz \
	    $gt \
	    chr$chr \
	    $binsize
        rm $dir/$type.$norm.chr$chr.txt
    done
#    for type in expected norm
#    do
#        $juicertool dump $type $norm $hic.hic $chr BP $binsize $dir/$type.$norm.chr$chr.matrix -d
#    done
done
