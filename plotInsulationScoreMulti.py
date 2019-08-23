#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from HiCmodule import JuicerMatrix
from InsulationScore import getInsulationScoreOfMultiSample
from generateCmap import *
from PlotModule import pltxticks

#import pdb

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    tp = lambda x:list(map(str, x.split(':')))
    parser.add_argument("input", help="<Input direcoty>:<label>", type=tp, nargs='*')
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="chromosome", type=str)
    parser.add_argument("--type", help="normalize type", type=str, default="KR")
    parser.add_argument("-r", "--resolution", help="resolution", type=int, default=25000)
    parser.add_argument("-s", "--start", help="start bp", type=int, default=0)
    parser.add_argument("-e", "--end", help="end bp", type=int, default=1000000)
    parser.add_argument("-x", "--sizex", help="Size of x axis (default: 30)", type=int, default=30)
    parser.add_argument("-y", "--sizey", help="Size of y axis (default: 2)", type=int, default=2)

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
        oe = dir + "/matrix/intrachromosomal/" + str(resolution) + "/oe."  + type + "." + chr + ".matrix.gz"
        eigen = dir + "/eigen/" + str(resolution) + "/gd_eigen."  + type + "." + chr + ".txt"

        print (observed)
        print (oe)
        print (eigen)

        samples.append(JuicerMatrix("RPM", observed, oe, eigen, resolution))
        print ("\n")

    Matrix = getInsulationScoreOfMultiSample(samples, labels)

#    pdb.set_trace()

    plt.figure(figsize=(11,2))
    plt.imshow(Matrix.T.iloc[:,s:e], cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']), aspect="auto")
    plt.colorbar()
    pltxticks(0, e-s, figstart, figend, 20)
    plt.yticks(np.arange(len(labels)), labels)

    plt.savefig(args.output + ".png")
