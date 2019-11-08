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

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx16384m -jar $pwd/../binaries/Aidenlab/juicer_tools.1.8.9_jcuda.0.8.jar"
gt=$($pwd/../script_rnakato/database.sh)/UCSC/$build/genome_table
chrlist=$($pwd/../script_rnakato/getchr_from_genometable.sh $gt)

dir=$matrixdir/Eigen/$binsize
mkdir -p $dir

ex(){ echo $1; eval $1; }

getPearson(){
    ex "$juicertool pearsons    -p $norm $hic chr$chr BP $binsize $dir/pearson.$norm.chr$chr.matrix"
    if test -s $dir/pearson.$norm.chr$chr.matrix; then
	gzip -f $dir/pearson.$norm.chr$chr.matrix
    fi
}

getEigen(){
    ex "$juicertool eigenvector -p $norm $hic chr$chr BP $binsize $dir/eigen.$norm.chr$chr.txt"
    $pwd/fixEigendir.py $dir/eigen.$norm.chr$chr.txt \
			$dir/eigen.$norm.chr$chr.txt.temp \
			$($pwd/../script_rnakato/database.sh)/UCSC/$build/refFlat.txt \
			chr$chr \
			$binsize
    mv $dir/eigen.$norm.chr$chr.txt.temp $dir/eigen.$norm.chr$chr.txt
    gzip -f $dir/eigen.$norm.chr$chr.txt
}

for chr in $chrlist
do
    if test $chr = "chrY" -o $chr = "chrM" -o $chr = "chrMT" ;then continue; fi

    chr=$(echo $chr | sed -e 's/chr//g')
    echo "chr$chr"
    if test ! -e $dir/pearson.$norm.chr$chr.matrix.gz; then
	getPearson &
    fi

    if test ! -e $dir/eigen.$norm.chr$chr.txt.gz; then
	getEigen &
    fi
done
