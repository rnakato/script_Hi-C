#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from HiCmodule import JuicerMatrix
from InsulationScore import getInsulationScoreOfMultiSample
from generateCmap import *
from loadData import *
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
    parser.add_argument("-r", "--resolution", help="resolution", type=int, default=25000)
    parser.add_argument("-s", "--start", help="start bp", type=int, default=0)
    parser.add_argument("-e", "--end", help="end bp", type=int, default=1000000)
    parser.add_argument("--vmax", help="max value of color bar", type=int, default=50)
    parser.add_argument("--vmin", help="min value of color bar", type=int, default=0)
    parser.add_argument("-d", "--vizdistancemax", help="max distance in heatmap", type=int, default=0)
#    parser.add_argument("--xsize", help="xsize for figure", type=int, default=10)
    parser.add_argument("--ysize", help="ysize (* times of samples)", type=int, default=3)

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
    if (args.log):
        vmax = np.log1p(vmax)
        vmin = np.log1p(vmin)
    figsize_x = max(int((figend-figstart)/2000000), 10)

    print (chr)
    print (resolution)
    samples = []
    for dir in dirs:
        observed = dir + "/Matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
        eigen = dir + "/eigen/" + str(resolution) + "/gd_eigen."  + type + "." + chr + ".txt"
        samples.append(JuicerMatrix("RPM", observed, resolution, eigenfile=eigen))

    ### Plot
    nsample = len(samples)
    plt.figure(figsize=(figsize_x, nsample*args.ysize))

    for i, sample in enumerate(samples):
        tadfile = dirs[i] + "/TAD/" + type + "/" + str(resolution) + "_blocks.bedpe"
        print(tadfile)
        tads = loadTADs(tadfile, chr[3:], start=figstart, end=figend)
        loopfile = dirs[i] + "/loops/" + type + "/merged_loops.bedpe"
        print(loopfile)
        loops = loadloops(loopfile, chr[3:], start=figstart, end=figend)
        # Hi-C Map
        plt.subplot2grid((nsample*2, 4), (i*2,0), rowspan=2, colspan=4)
        if (args.log):
            drawHeatmapTriangle(plt, sample.getlog(), resolution,
                                figstart=figstart, figend=figend,
                                tads=tads, loops=loops,
                                vmax=vmax, vmin=vmin, distancemax=args.vizdistancemax,
                                label=labels[i], xticks=True)
        else:
            drawHeatmapTriangle(plt, sample.getmatrix(), resolution,
                                figstart=figstart, figend=figend,
                                tads=tads, loops=loops,
                                vmax=vmax, vmin=vmin, distancemax=args.vizdistancemax,
                                label=labels[i], xticks=True)

    plt.tight_layout()
    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
