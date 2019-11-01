Aidendir=/work/git/binaries/Aidenlab/
Datadir=/work/Database
ncore=32

jdir=$Aidendir/juicer_1.5.7  ####  <= new one!!!
juicertool="java -Xmx16384m -jar $Aidendir/juicer_tools.1.8.9_jcuda.0.8.jar"

odir=$1
label=$(basename $odir)
build=$2
enzyme=$3
fastq_post=$4
gt=$Datadir/UCSC/$build/genome_table
bwaindex=$Datadir/bwa-indexes/UCSC-$build

bash $jdir/CPU/juicer.sh -t $ncore -g $build -d $odir \
     -s $enzyme -a $label -p $gt \
     -z $bwaindex -D $jdir -e $fastq_post -S map

bash $jdir/CPU/juicer.sh -t $ncore -g $build -d $odir \
     -s $enzyme -a $label -p $gt \
     -z $bwaindex -D $jdir -e $fastq_post -S aftermap
