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
juicertool="java -Xms512m -Xmx2048m -jar $pwd/../binaries/Aidenlab/juicer_tools.1.8.9_jcuda.0.8.jar"
gt=$($pwd/../script_rnakato/database.sh)/UCSC/$build/genome_table
chrlist=$($pwd/../script_rnakato/getchr_from_genometable.sh $gt)

dir=$matrixdir/Eigen/$binsize
mkdir -p $dir

ex(){ echo $1; eval $1; }

for chr in $chrlist
do
    if test $chr = "chrY" -o $chr = "chrM" -o $chr = "chrMT" ;then continue; fi

    chr=$(echo $chr | sed -e 's/chr//g')
    echo "chr$chr"
    ex "$juicertool pearsons    -p $norm $hic chr$chr BP $binsize $dir/pearson.$norm.chr$chr.matrix"
    ex "$juicertool eigenvector -p $norm $hic chr$chr BP $binsize $dir/eigen.$norm.chr$chr.txt"
    if test -s $dir/pearson.$norm.chr$chr.matrix; then
	gzip -f $dir/pearson.$norm.chr$chr.matrix
    fi
    if test -s $dir/eigen.$norm.chr$chr.txt; then
	$pwd/fixEigendir.py $dir/eigen.$norm.chr$chr.txt \
			    $dir/eigen.$norm.chr$chr.txt.temp \
			    $($pwd/../script_rnakato/database.sh)/UCSC/$build/refFlat.txt \
			    chr$chr \
			    $binsize
	mv $dir/eigen.$norm.chr$chr.txt.temp $dir/eigen.$norm.chr$chr.txt
	gzip -f $dir/eigen.$norm.chr$chr.txt
    fi
done
