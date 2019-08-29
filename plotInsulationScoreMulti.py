#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import ndimage
from HiCmodule import JuicerMatrix
from InsulationScore import getInsulationScoreOfMultiSample
from generateCmap import *
from PlotModule import pltxticks

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
#        eigen = ""
 #       print (observed)
 #       print (oe)
#        print (eigen)

        samples.append(JuicerMatrix("RPM", observed, eigen, resolution))
    #        print ("\n")

    nrow_heatmap = 3
    nrow_eigen = 1
    nrow_now = 0

    ### Plot
    if (args.multi):
        plt.figure(figsize=(10, 14))
        nrow = nrow_heatmap + nrow_eigen + len(samples)
    else:
        plt.figure(figsize=(10,6))
        nrow = nrow_heatmap + nrow_eigen + 2

    # Hi-C Map
    plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=nrow_heatmap, colspan=4)
    dst = ndimage.rotate(samples[0].getmatrix().iloc[s:e,s:e], 45,
                         order=0, reshape=True, prefilter=False, cval=0)
    img = plt.imshow(dst, clim=(0, 50), cmap=generate_cmap(['#FFFFFF', '#d10a3f']),
                     interpolation="nearest", aspect='auto')
    plt.ylim(int(dst.shape[0]/2)+1,0)
    plt.title(labels[0])
    #        pltxticks(0, (e-s)*1.41, figstart, figend, 10)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False  # labels along the bottom edge are off
    )

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
        plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=2, colspan=5)
        Matrix = getInsulationScoreOfMultiSample(samples, labels)
        plt.imshow(Matrix.T.iloc[:,s:e],
                   clim=(0.4, 1.0),
                   cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']),
                   aspect="auto")
        plt.colorbar()
        pltxticks(0, e-s, figstart, figend, 10)
        plt.yticks(np.arange(len(labels)), labels)

    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
