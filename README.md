# script_Hi-C
scripts for Hi-C analysis

#### distance_vs_count.Juicer
Count genomic distance of read pairs in the input file (supposing align/merged_nodups.txt.gz in Juicer output.)

Usage: 

      distance_vs_count.Juicer <file> <winsize> <MAPQ>
            <file>:    Input file  (merged_nodups.txt.gz)
            <winsize>: window size (default: 10000)
            <MAPQ>:    MAPQ threshold (default: 30)

#### convert_JuicerDump_to_dense.py
Convert interaction frequency file dumped by Juicer to dense matrix

Usage: 

      python /home/git/script_Hi-C/convert_JuicerDump_to_dense.py <inputfile> <outputfile> <genometable> <chr> <resolution> [--help]


#### plotHiCMatrix.py
Plot Hi-C intercation heatmap from Juicer matrix

Usage:

     plotHiCMatrix.py <matrix> <output name (png)> <start> <end> <title in figure>
     Example: plotHiCMatrix.py WT/intrachromosomal/25000/observed.KR.chr7.matrix.gz chr7/WT_chr7_25000000-31000000.png 25000000 31000000 WT

                                 
#### plotInsulationScore.py
Plot insulation score from Juicer matrix

Usage:

     plotInsulationScore.py [-h] [--num4norm NUM4NORM] [--distance DISTANCE]
                                 [--sizex SIZEX] [--sizey SIZEY]
                                 matrix output resolution

#### plotMultiScaleInsulationScore.py
Plot multi-scale insulation scores from Juicer matrix

Usage:

     plotMultiScaleInsulationScore.py [-h] [--num4norm NUM4NORM]
                                           [--sizex SIZEX] [--sizey SIZEY]
                                           matrix output resolution

