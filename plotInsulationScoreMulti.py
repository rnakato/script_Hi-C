#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from HiCmodule import JuicerMatrix
from InsulationScore import getInsulationScoreOfMultiSample
from generateCmap import *

import pdb


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Direcotyr of input matrix", type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="chromosome", type=str)
    parser.add_argument("--type", help="normalize type", type=str, default="KR")
    parser.add_argument("--resolution", help="resolution", type=int, default=25000)
    parser.add_argument("--start", help="start bp", type=int, default=0)
    parser.add_argument("--end", help="end bp", type=int, default=1000000)
    parser.add_argument("--sizex", help="Size of x axis (default: 30)", type=int, default=30)
    parser.add_argument("--sizey", help="Size of y axis (default: 2)", type=int, default=2)

    args = parser.parse_args()
    #    print(args)

    dir = args.directory
    chr = args.chr
    resolution = args.resolution
    type = args.type
    figstart = args.start
    figend = args.end
    s = int(figstart / resolution)
    e = int(figend   / resolution)

    observed = dir + "/matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
    oe = dir + "/matrix/intrachromosomal/" + str(resolution) + "/oe."  + type + "." + chr + ".matrix.gz"
    eigen = dir + "/eigen/" + str(resolution) + "/gd_eigen."  + type + "." + chr + ".txt"

    print (observed)
    print (oe)
    print (eigen)
    print (chr)
    print (resolution)

    print("a = JuicerMatrix('RPM', observed, oe, eigen, resolution)")
    a = JuicerMatrix("RPM", observed, oe, eigen, resolution)

    Matrix = getInsulationScoreOfMultiSample([a, a, a], ["a", "b", "c"])

#    pdb.set_trace()


    plt.figure(figsize=(11,2))
    plt.imshow(Matrix.T.iloc[:,s:e], cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']), aspect="auto")
    plt.colorbar()
    plt.savefig(args.output + ".png")
#    pltxticks(0, e-s, figstart, figend, 20)
 #   plt.yticks(np.arange(nsample), labels_full)

    # generate MI .png
#    fig, ax = plt.subplots(1, 1, figsize=(args.sizex, args.sizey))
 #   plt.imshow(MI.MI, clim=(-1, 1), cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']), aspect="auto")
  #  plt.colorbar()
   # plt.savefig(args.output + ".png")
