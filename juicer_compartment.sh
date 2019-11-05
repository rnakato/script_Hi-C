#!/bin/bash

norm=$1
hic=$2
odir=$3
build=$4

for res in 25000 50000 100000; do
    makeEigen.sh $norm $odir $hic $res $build
done
