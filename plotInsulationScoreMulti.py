#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from HiCmodule import JuicerMatrix
from InsulationScore import getInsulationScoreOfMultiSample
from generateCmap import *
from PlotModule import *
from loadData import *

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
    parser.add_argument("--multi", help="plot MultiInsulation Score", action='store_true')
    parser.add_argument("--vmax", help="max value of color bar", type=int, default=50)
    parser.add_argument("--vmin", help="min value of color bar", type=int, default=0)

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
        observed = dir + "/matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
        eigen = dir + "/eigen/" + str(resolution) + "/gd_eigen."  + type + "." + chr + ".txt"
 #       print (observed)
#        print (eigen)

        samples.append(JuicerMatrix("RPM", observed, resolution, eigenfile=eigen))

    nrow_heatmap = 3
    nrow_eigen = 1
    nrow_now = 0

    ### Plot
    if (args.multi):
        plt.figure(figsize=(10, 14))
        nrow = nrow_heatmap + nrow_eigen + len(samples)
    else:
        plt.figure(figsize=(10,6))
        nrow = nrow_heatmap + nrow_eigen + 4

    # Hi-C Map
    plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=nrow_heatmap, colspan=5)

    tadfile = dirs[0] + "/contact_domain/" + str(resolution) + "_blocks.bedpe"
    print(tadfile)
    tads = loadTADs(tadfile, chr[3:], start=figstart, end=figend)
    loopfile = dirs[0] + "/loops/merged_loops.bedpe"
    print(loopfile)
    loops = loadloops(loopfile, chr[3:], start=figstart, end=figend)

    drawHeatmapTriangle(plt, samples[0].getmatrix(), resolution,
                        figstart=figstart, figend=figend,
                        tads=tads, loops=loops,
                        vmax=vmax, vmin=vmin,
                        label=labels[0], xticks=False)

    nrow_now += nrow_heatmap

    # Compartment
    plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=nrow_eigen, colspan=4)
    plt.plot(samples[0].getEigen())
    plt.xlim([s,e])
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False  # labels along the bottom edge are off
    )

    nrow_now += nrow_eigen

    # Insulation score
    if (args.multi):
        for i, sample in enumerate(samples):
            plt.subplot2grid((nrow, 5), (i + nrow_now, 0), rowspan=1, colspan=5)
            MI = sample.getMultiInsulationScore()
            plt.imshow(MI.T.iloc[:,s:e],
                       clim=(0.4, 1.0),
                       cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']),
                       aspect="auto")
            plt.title(labels[i])
            pltxticks(0, e-s, figstart, figend, 10)
#            plt.yticks([2], (labels[i]))

#            plt.xlim([s,e])
#            pltxticks(s, e, figstart, figend, 10)
            plt.colorbar()
        plt.tight_layout()
    else:
        Matrix = getInsulationScoreOfMultiSample(samples, labels)
        plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=2, colspan=5)
        plt.imshow(Matrix.T.iloc[:,s:e],
                   clim=(0.4, 1.0),
                   cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']),
                   aspect="auto")
        plt.colorbar()
        pltxticks(0, e-s, figstart, figend, 10)
        plt.yticks(np.arange(len(labels)), labels)

        nrow_now += 2

        plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=2, colspan=4)
        for i in range(len(samples)):
            plt.plot(Matrix.iloc[s:e,i], label=labels[i])
        plt.xlim([figstart, figend])
#        plt.legend()

    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
