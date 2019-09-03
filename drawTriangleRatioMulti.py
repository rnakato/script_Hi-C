#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import ndimage
from HiCmodule import JuicerMatrix
from generateCmap import *
from Cluster import make3dmatrixRatio
from PlotModule import *

def getDirectionalFreqRatio(mat, resolution, strand, *, distance=2000000):
    arraysize = mat.shape[0]
    array = np.zeros(arraysize)
    nbin = int(distance/resolution)
    for i in range(nbin, arraysize - nbin):
        if (strand == "+"):
            val = mat[i+1:i+nbin, i].mean()
        else:
            val = mat[i, i-nbin:i-1].mean()
        array[i] = val

    return array

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
    parser.add_argument("--vmax", help="max value of color bar", type=int, default=2)
    parser.add_argument("--vmin", help="min value of color bar", type=int, default=-2)

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

    ### Plot
    smooth_median_filter = 3
    EnrichMatrices = make3dmatrixRatio(samples, smooth_median_filter)

    nsample = len(samples) -1
    plt.figure(figsize=(10, nsample*5))

    for i, sample in enumerate(EnrichMatrices):
        # Hi-C Map
        plt.subplot2grid((nsample*4, 5), (i*4,0), rowspan=2, colspan=5)
        drawHeatmapTriangle(plt, sample, resolution,
                            figstart=figstart, figend=figend, vmax=vmax, vmin=vmin,
                            cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
                            distancemax=args.vizdistancemax,
                            label=labels[i+1], xticks=True)

        DFRplus = getDirectionalFreqRatio(sample, resolution, "+")
        DFRminus = getDirectionalFreqRatio(sample, resolution, "-")
        plt.subplot2grid((nsample*4, 5), (i*4+2,0), rowspan=1, colspan=4)
        plt.plot(DFRplus, label="Right")
        plt.plot(DFRminus, label="Left")
        plt.xlim([s,e])
        plt.legend()
        plt.subplot2grid((nsample*4, 5), (i*4+3,0), rowspan=1, colspan=4)
        plt.bar(range(len(DFRplus)), DFRplus - DFRminus)
        plt.xlim([s,e])

    plt.tight_layout()
    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
