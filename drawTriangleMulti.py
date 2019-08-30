#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from HiCmodule import JuicerMatrix
from InsulationScore import getInsulationScoreOfMultiSample
from generateCmap import *
from loadData import loadTADs
from PlotModule import *

#import pdb

def main():
    parser = argparse.ArgumentParser()
    tp = lambda x:list(map(str, x.split(':')))
    parser.add_argument("input", help="<Input direcoty>:<label>", type=tp, nargs='*')
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="chromosome", type=str)
    parser.add_argument("--log", help="logged count", action='store_true')
    parser.add_argument("--type", help="normalize type", type=str, default="KR")
    parser.add_argument("-d", "--distance", help="distance for DI", type=int, default=500000)
    parser.add_argument("-r", "--resolution", help="resolution", type=int, default=25000)
    parser.add_argument("-s", "--start", help="start bp", type=int, default=0)
    parser.add_argument("-e", "--end", help="end bp", type=int, default=1000000)

    args = parser.parse_args()
#    print(args)

    dirs = []
    labels = []
    for input in args.input:
        dirs.append(input[0])
        if (len(input) >1):
            labels.append(input[1])
        else:
            labels.append("")

    chr = args.chr
    resolution = args.resolution
    type = args.type
    figstart = args.start
    figend = args.end
    s = int(figstart / resolution)
    e = int(figend   / resolution)

    print (chr)
    print (resolution)
    samples = []
    for dir in dirs:
        observed = dir + "/matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
#        oe = dir + "/matrix/intrachromosomal/" + str(resolution) + "/oe."  + type + "." + chr + ".matrix.gz"
        eigen = dir + "/eigen/" + str(resolution) + "/gd_eigen."  + type + "." + chr + ".txt"
        #        print (observed)
        #print (eigen)

        samples.append(JuicerMatrix("RPM", observed, resolution, eigenfile=eigen))

    ### Plot

    nsample = len(samples)
    plt.figure(figsize=(10, nsample*3))

    for i, sample in enumerate(samples):
        tadfile = dirs[i] + "/contact_domain/" + str(resolution) + "_blocks.bedpe"
        print(tadfile)
        tads = loadTADs(tadfile, chr[3:], start=figstart, end=figend)
        # Hi-C Map
        plt.subplot2grid((nsample*2, 4), (i*2,0), rowspan=2, colspan=4)
        if (args.log):
            drawHeatmapTriangle(plt, sample.getlog(), resolution,
                                figstart=figstart, figend=figend, tads=tads,
                                vmax=6, label=labels[i], xticks=False)
        else:
            drawHeatmapTriangle(plt, sample.getmatrix(), resolution,
                                figstart=figstart, figend=figend, tads=tads,
                                vmax=50, label=labels[i], xticks=False)
#        dst = ndimage.rotate(df.iloc[s:e,s:e], 45,
#                             order=0, reshape=True, prefilter=False, cval=0)
#        img = plt.imshow(dst, clim=(0, valmax), cmap=generate_cmap(['#FFFFFF', '#d10a3f']),
#                         interpolation="nearest", aspect='auto')
#        plt.ylim(int(dst.shape[0]/2)+1,0)
#        plt.title(labels[i])
        #        pltxticks(0, (e-s)*1.41, figstart, figend, 10)

    plt.tight_layout()
    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
