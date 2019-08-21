#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <norm> <outputdir> <hic file>" 1>&2
}

if [ $# -ne 3 ]; then
  usage
  exit 1
fi

norm=$1
odir=$2
hic=$3

ex(){ echo $1; eval $1; }

pwd=$(cd $(dirname $0) && pwd)
juicertool="java -Xms512m -Xmx2048m -jar $pwd/../binaries/Aidenlab/juicer_tools.1.8.9_jcuda.0.8.jar"

hicdir=$odir/loops/$norm
mkdir -p $hicdir
ex "$juicertool hiccups -r 5000,10000,25000 -k $norm --ignore_sparsity $hic $hicdir"

# APA
file=$odir/loops/$norm/merged_loops.bedpe
if test -s $file; then
    mkdir -p $odir/APA
    ex "$juicertool apa $hic $file $odir/APA"
else
    echo "APA Warning: $file does not exist."
fi
