# script_Hi-C
scripts for Hi-C analysis

#### distance_vs_count.Juicer
Count genomic distance of read pairs in the input file (supposing align/merged_nodups.txt.gz in Juicer output.)

Usage: 

      distance_vs_count.Juicer <file> <winsize> <MAPQ>
            <file>:    Input file  (merged_nodups.txt.gz)
            <winsize>: window size (default: 10000)
            <MAPQ>:    MAPQ threshold (default: 30)

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

