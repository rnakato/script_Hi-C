#!/bin/bash

norm=$1
odir=$2
build=$3

for res in 100000 50000 25000
do
    makeInslationScore.sh $norm $odir $res $build
done
