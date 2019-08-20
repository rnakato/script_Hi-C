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

if test $build = "mm10" -o $build = "mm9"; then
    chrnum=19
else
    chrnum=22
fi

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx2048m -jar $pwd/../binaries/Aidenlab/juicer_tools.1.8.9_jcuda.0.8.jar"
dir=$matrixdir/Eigen/$binsize
mkdir -p $dir

pwd=$(cd $(dirname $0) && pwd)

ex(){
    command=$1
    echo $1
    `echo $1`
}

for chr in $(seq 1 $chrnum) X; do
    echo $chr
    ex "$juicertool pearsons    -p $norm $hic chr$chr BP $binsize $dir/pearson.$norm.chr$chr.matrix"
    ex "$juicertool eigenvector -p $norm $hic chr$chr BP $binsize $dir/eigen.$norm.chr$chr.txt"
    gzip -f $dir/pearson.$norm.chr$chr.matrix $dir/eigen.$norm.chr$chr.txt
done
