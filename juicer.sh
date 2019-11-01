Datadir=/work/Database
ncore=32

build=$2
enzyme=$3
fastq_post=$4


ex(){ echo $1; eval $1; }

doBWA(){
    label=$1
    bash $jdir/CPU/juicer.sh -t $ncore -g $build -d $odir -s $enzyme -a $label -p $gt \
	 -z $bwaindex -D $jdir -e $fastq_post -S map
}

doPostProcessing(){
    label=$1
    bash $jdir/CPU/juicer.sh -t $ncore -g $build -d $odir -s $enzyme -a $label -p $gt \
	 -z $bwaindex -D $jdir -e $fastq_post -S aftermap
}

dopigz(){
    pigz $odir/*/*.sam $odir/splits/SRR5266587.fastq.gz.sort.txt
    pigz $odir/aligned/merged_nodups.txt $odir/aligned/merged_sort.txt
}

dounpigz(){
    unpigz $odir/*/*.sam.gz $odir/splits/SRR5266587.fastq.gz.sort.txt.gz
    unpigz $odir/aligned/merged_nodups.txt.gz $odir/aligned/merged_sort.txt.gz
}

call_TAD(){
    norm=$1
    hic=$2

    dir=$odir/TAD/$norm
    mkdir -p $dir
    for res in 10000 25000 50000; do
	ex "$juicertool arrowhead --ignore_sparsity -m 2000 -r $res -k $norm $hic $dir"
#	cat $dir/${res}_blocks.bedpe | awk '{OFS="\t"} NR>1 {print $1, $2, $3 }' | sort -k1 -k2,2n > tmp1.txt
#	perl /work/sakata/Perl/border.pl tmp1.txt $res | sort -k1 -k2,2n | uniq > tmp2.txt
#	addchr.pl tmp2.txt > $dir/borders_${res}_blocks
#	addchr.pl $dir/${res}_blocks.bedpe | awk '{OFS="\t"} NR>1 {print }' > $dir/${res}_blocks.txt
#	rm -f tmp*.txt
    done
}

findmotif(){
    ex "$juicertool motifs $build $odir/loops/$norm"
}

make_matrix(){
    norm=$1
    hic=$2
    for res in 25000 50000 100000
    do
	ex "makeMatrix_intra.sh $norm $odir $hic $res $build"
	#    ex "makeMatrix_inter.sh $matrixdir "$hic" $res $build 0.0"
    done
}

call_compartment(){
    norm=$1
    hic=$2
    for res in 25000 50000 100000
    do
	ex "makeEigen.sh $norm $odir $hic $res $build"
    done
}

call_InsulationScore(){
    norm=$1
    for res in 100000 50000 25000
    do
	ex "makeInslationScore.sh $norm $odir $res $build"
    done
}

PlotDistance(){
    mkdir -p $odir/distance
    ex "distance_vs_count.Juicer $odir/aligned/merged_nodups.txt.gz > $odir/distance/distance_vs_count.10kb.MAPQ30.txt"
}
