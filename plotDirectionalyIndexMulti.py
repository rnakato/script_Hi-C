#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
#from scipy import ndimage
from HiCmodule import JuicerMatrix
from DirectionalityIndex import getDirectionalityIndexOfMultiSample
from generateCmap import *
from loadData import loadTADs
from PlotModule import *
#import pysnooper

#import pdb

#@pysnooper.snoop()
def main():
    parser = argparse.ArgumentParser()
    tp = lambda x:list(map(str, x.split(':')))
    parser.add_argument("input", help="<Input direcoty>:<label>", type=tp, nargs='*')
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="chromosome", type=str)
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
    length = figend - figstart
    binnum = e-s

    print ("width: " + str(length) + ", " + str(binnum) + " bins.")
    if (length <= 0):
        print ("Error: end < start.")
        exit(1)
    print (chr)
    print (resolution)

    samples = []
    for dir in dirs:
        observed = dir + "/matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
        eigen = dir + "/eigen/" + str(resolution) + "/gd_eigen."  + type + "." + chr + ".txt"
#        print (observed)
 #       print (eigen)
        samples.append(JuicerMatrix("RPM", observed, resolution, eigenfile=eigen))

    ### Plot
    plt.figure(figsize=(10,6))

    # Hi-C Map
    plt.subplot2grid((7, 5), (0,0), rowspan=3, colspan=4)
    tadfile = dirs[0] + "/contact_domain/" + str(resolution) + "_blocks.bedpe"
    print(tadfile)
    tads = loadTADs(tadfile, chr[3:], start=figstart, end=figend)

    drawHeatmapTriangle(plt, samples[0].getmatrix(), resolution,
                        figstart=figstart, figend=figend, tads=tads,
                        vmax=50, label=labels[0], xticks=False)

    # Compartment
    plt.subplot2grid((7, 5), (3,0), rowspan=1, colspan=4)
    plt.plot(samples[0].getEigen())
    plt.xlim([s,e])
    xtickoff(plt)

    # DI
    plt.subplot2grid((7, 5), (4,0), rowspan=3, colspan=5)
    vDI = getDirectionalityIndexOfMultiSample(samples, labels, distance=args.distance)
    plt.imshow(vDI.iloc[:,s:e],
               clim=(-1000, 1000),
               cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
               aspect="auto")
    plt.colorbar()
    pltxticks(0, e-s, figstart, figend, 10)
    plt.yticks(np.arange(len(labels)), labels)

    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
