#!/bin/bash

ex(){ echo $1; eval $1; }

odir=$1
ex "unpigz $odir/*/*.sam.gz $odir/splits/*.fastq.gz.sort.txt.gz $odir/aligned/merged_nodups.txt.gz $odir/aligned/merged_sort.txt.gz"
