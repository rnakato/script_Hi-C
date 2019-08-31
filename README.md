# script_Hi-C
scripts for Hi-C analysis

## QC
### distance_vs_count.Juicer
Count genomic distance of read pairs in the input file (supposing align/merged_nodups.txt.gz in Juicer output.)

Usage: 

      distance_vs_count.Juicer <file> <winsize> <MAPQ>
            <file>:    Input file  (merged_nodups.txt.gz)
            <winsize>: window size (default: 10000)
            <MAPQ>:    MAPQ threshold (default: 30)

## Data generation
### convert_JuicerDump_to_dense.py
Convert interaction frequency file dumped by Juicer to dense matrix

Usage: 

    convert_JuicerDump_to_dense.py <inputfile> <outputfile> <genometable> <chr> <resolution> [--help]

Example:

    convert_JuicerDump_to_dense.py \
        Matrix.observed.VC_SQRT.chrX.txt \
        Matrix.observed.VC_SQRT.chrX.matrix.gz \
        genome_table.txt \
        chrX \
        25000


### makeMatrix_intra.sh
Generate dense matrix from .hic for each chromosome

Usage:

     makeMatrix_intra.sh <normalize type (e.g. KR)> <output directory> <.hic> <resolution> <build (r.g., hg38)>

### makeEigen.sh
Generate eigenvector file in that +- of the value is adjusted by the number of genes

Usage:

     makeEigen.sh <normalize type (e.g. KR)> <output directory> <.hic> <resolution> <build (r.g., hg38)>

### makeInslationScore.sh
Calculate insulation score file from dense matrix data

Usage:

     makeInslationScore.sh <normalize type (e.g. KR)> <output directory> <resolution> <build (r.g., hg38)>

## Visualization
### plotHiCMatrix.py
Plot Hi-C intercation heatmap from Juicer matrix

Usage:

    plotHiCMatrix.py <matrix> <output name (png)> <start> <end> <title in figure>

Example: 

    plotHiCMatrix.py WT/intrachromosomal/25000/observed.KR.chr7.matrix.gz chr7/WT_chr7_25000000-31000000.png 25000000 31000000 WT
                                 
### plotInsulationScore.py
Plot insulation score from Juicer matrix

Usage:

     plotInsulationScore.py [-h] [--num4norm NUM4NORM] [--distance DISTANCE]
                                 [--sizex SIZEX] [--sizey SIZEY]
                                 matrix output resolution

### plotMultiScaleInsulationScore.py
Plot multi-scale insulation scores from Juicer matrix

Usage:

     plotMultiScaleInsulationScore.py [-h] [--num4norm NUM4NORM]
                                           [--sizex SIZEX] [--sizey SIZEY]
                                           matrix output resolution


### plotDirectionalyIndexMulti.py
Plot multi-scale insulation scores from Juicer matrix

Usage:

     plotDirectionalyIndexMulti.py $samples $head.DI $region --type $type
### plotInsulationScoreMulti.py
Plot multi-scale insulation scores from Juicer matrix

Usage:

     plotInsulationScoreMulti.py $samples $head.MultiIS $region --type $type --multi
     plotInsulationScoreMulti.py $samples $head.singleIS $region --type $type
### drawTriangleMulti.py
Plot multi-scale insulation scores from Juicer matrix

Usage:

     drawTriangleMulti.py $samples $head.Heatmap $region --type $type
     drawTriangleMulti.py $samples $head.Heatmaplog $region --type $type --log
### drawTriangleRatioMulti.py 
Plot multi-scale insulation scores from Juicer matrix

Usage:

     drawTriangleRatioMulti.py $samples $head.HeatmapRatio $region --type $type
