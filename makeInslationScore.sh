#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <norm> <matrixdir> <binsize> <build>" 1>&2
}

if [ $# -ne 4 ]; then
  usage
  exit 1
fi

norm=$1
matrixdir=$2
binsize=$3
build=$4

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx2048m -jar $pwd/../binaries/Aidenlab/juicer_tools.1.8.9_jcuda.0.8.jar"
gt=$($pwd/../script_rnakato/database.sh)/UCSC/$build/genome_table
chrlist=$($pwd/../script_rnakato/getchr_from_genometable.sh $gt)

dir=$matrixdir/InsulationScore/$binsize
mkdir -p $dir

ex(){ echo $1; eval $1; }

for chr in $chrlist
do
    if test $chr = "chrY" -o $chr = "chrM" -o $chr = "chrMT" ;then continue; fi

    i=$(echo $chr | sed -e 's/chr//g')
    echo "chr$i"
    matrix=$matrixdir/Matrix/intrachromosomal/$binsize/observed.$norm.chr$i.matrix.gz

    if test -s $matrix; then
	InsulationScore.py               $matrix $dir/chr$i chr$i $resolution
	plotInsulationScore.py           $matrix $dir/chr$i       $resolution
	plotMultiScaleInsulationScore.py $matrix $dir/chr$i       $resolution
    else
	echo "Warning: $matrix does not exist."
    fi
done
