#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import ndimage
from HiCmodule import JuicerMatrix
from generateCmap import *
from PlotModule import *
from DirectionalFreqRatio import *

#import pdb
def main():
    parser = argparse.ArgumentParser()
    tp = lambda x:list(map(str, x.split(':')))
    parser.add_argument("input", help="<Input direcoty>:<label>", type=tp, nargs='*')
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="chromosome", type=str)
    parser.add_argument("--type", help="normalize type", type=str, default="KR")
    parser.add_argument("-r", "--resolution", help="resolution", type=int, default=25000)
    parser.add_argument("-s", "--start", help="start bp", type=int, default=0)
    parser.add_argument("-e", "--end", help="end bp", type=int, default=1000000)
    parser.add_argument("-d", "--vizdistancemax", help="max distance in heatmap", type=int, default=0)
    parser.add_argument("--xsize", help="xsize for figure", type=int, default=10)
 #   parser.add_argument("--ysize", help="ysize (* times of samples)", type=int, default=5)
    parser.add_argument("--vmax", help="max value of color bar", type=int, default=1)
    parser.add_argument("--vmin", help="min value of color bar", type=int, default=-1)

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
    vmax = args.vmax
    vmin = args.vmin

    print (chr)
    print (resolution)
    samples = []
    for dir in dirs:
        observed = dir + "/Matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
        samples.append(JuicerMatrix("RPM", observed, resolution))

    smooth_median_filter = 3
    EnrichMatrices = make3dmatrixRatio(samples, smooth_median_filter)

    rowspan_heatmap = 3
    rowspan_barplot = 1
    rowspan = rowspan_heatmap + rowspan_barplot*2

    nsample = len(samples) -1
    plt.figure(figsize=(args.xsize, nsample*rowspan))

    for i, sample in enumerate(EnrichMatrices):
        # Hi-C Map
        plt.subplot2grid((nsample*rowspan, 5), (i*rowspan,0), rowspan=rowspan_heatmap, colspan=5)
        drawHeatmapTriangle(plt, sample, resolution,
                            figstart=figstart, figend=figend, vmax=vmax, vmin=vmin,
                            cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
                            distancemax=args.vizdistancemax,
                            label=labels[i+1], xticks=True)

        dfr = DirectionalFreqRatio(sample, resolution)

        plt.subplot2grid((nsample*rowspan, 5), (i*rowspan + rowspan_heatmap,0), rowspan=rowspan_barplot, colspan=4)
        plt.plot(dfr.getarrayplus(), label="Right")
        plt.plot(dfr.getarrayminus(), label="Left")
        plt.xlim([s,e])
        plt.legend(loc='lower right')

        plt.subplot2grid((nsample*rowspan, 5), (i*rowspan + rowspan_heatmap + rowspan_barplot,0), rowspan=rowspan_barplot, colspan=4)
        diff = dfr.getarraydiff()
        plt.bar(range(len(diff)), diff)
        plt.xlim([s,e])
        plt.title("Right - Left")

    plt.tight_layout()
    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
