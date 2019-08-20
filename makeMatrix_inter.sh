#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <norm> <matrixdir> <hic file> <binsize> <build> <lim_pzero>" 1>&2
}

if [ $# -ne 5 ]; then
  usage
  exit 1
fi

norm=$1
dir=$2
hic=$3
binsize=$4
build=$5
lim_pzero=$6

gt=$(../script_rnakato/database.sh)/UCSC/$build/genome_table

if test $build = "mm10" -o $build = "mm9"; then
    chrnum=19
else
    chrnum=22
fi

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx2048m -jar $pwd/../binaries/Aidenlab/juicer_tools.1.8.9_jcuda.0.8.jar"

for i in $(seq 1 $chrnum)
do
    for j in $(seq $i $chrnum); do
        d=$dir/interchromosomal/$binsize/chr$i-chr$j
        mkdir -p $d
        for type in observed oe; do
            echo $i $j $type
            $juicertool dump $type VC_SQRT $hic $i $j BP $binsize $d/$type.matrix -d
	    gzip -f $d/$type.matrix
        done
    done
done

for str in observed #oe
do
    $pwd/merge_JuicerMatrix_to_Genome.py $dir/interchromosomal \
        $dir/interchromosomal/$binsize/genome.$str.full.$lim_pzero.pickle \
        $binsize $str $lim_pzero $chrnum
    $pwd/merge_JuicerMatrix_to_Genome.py $dir/interchromosomal \
        $dir/interchromosomal/$binsize/genome.$str.evenodd.$lim_pzero.pickle \
        $binsize $str $lim_pzero $chrnum --evenodd
done
