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

findmotif(){
    ex "$juicertool motifs $build $odir/loops/$norm"
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
