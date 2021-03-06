import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import ndimage
from generateCmap import *

Mega = 1000000

def pltxticks(start, end, figstart, figend, nmem):
    mem = int((end - start)/nmem)
    x = start + np.arange(nmem+1) * mem
    val = (figend - figstart) / nmem
    xval = (figstart + np.arange(nmem+1)*val)/Mega
    plt.xticks(x, xval)

def pltyticks(start, end, figstart, figend, nmem):
    mem = int((end - start)/nmem)
    x = start + np.arange(nmem+1) * mem
    val = (figend - figstart) / nmem
    xval = (figstart + np.arange(nmem+1)*val)/Mega
    plt.yticks(x, xval)

def axxticks(ax, start, end, figstart, figend, nmem):
    mem = int((end - start)/nmem)
    x = start + np.arange(nmem+1) * mem
    val = (figend - figstart) / nmem
    xval = (figstart + np.arange(nmem+1)*val)/Mega
    ax.set_xticks(x)
    ax.set_xticklabels(xval)

def ytickoff(plt):
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,        # ticks along the bottom edge are off
        right=False,       # ticks along the top edge are off
        labelleft=False
    )

def xtickoff(plt):
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False  # labels along the bottom edge are off
    )

def drawHeatmapSquare(plt, matrix, resolution, *, tads="", loops="",
                      figstart=0, figend=0, distancemax=0,
                      vmin=0, vmax=50, label="", xticks=True,
                      cmap=generate_cmap(['#FFFFFF', '#d10a3f'])):
    s = int(figstart / resolution)
    e = int(figend   / resolution)
    if (e==0): e = matrix.shape[0]

    if (isinstance(matrix, pd.DataFrame)):
        mat = matrix.iloc[s:e,s:e]
    else:
        mat = matrix[s:e,s:e]

    plt.imshow(mat, clim=(vmin, vmax), cmap=cmap, interpolation="nearest", aspect=1)

    if (isinstance(tads, pd.DataFrame)):
        for tad in tads.itertuples(name=None):
            x1 = int(tad[2])/resolution - s
            x2 = int(tad[3])/resolution - s
            plt.vlines(x1, 0, matrix.shape[0], linestyle='dashed', linewidth=0.3)
            plt.vlines(x2, 0, matrix.shape[0], linestyle='dashed', linewidth=0.3)

    if (label != ""): plt.title(label)
    if (xticks):
        pltxticks(0, e-s, figstart, figend, 10)
    else:
        xtickoff(plt)
        ytickoff(plt)

def drawHeatmapTriangle(plt, matrix, resolution, *, tads="", loops="",
                        figstart=0, figend=0, distancemax=5000000,
                        vmin=0, vmax=50, label="", xticks=True,
                        cmap=generate_cmap(['#FFFFFF', '#d10a3f'])):
    s = int(figstart / resolution)
    e = int(figend   / resolution)
    if (e==0): e = matrix.shape[0]

    if (isinstance(matrix, pd.DataFrame)):
        mat = matrix.iloc[s:e,s:e]
    else:
        mat = matrix[s:e,s:e]

    dst = ndimage.rotate(mat, 45, order=0, reshape=True, prefilter=False, cval=0)
    plt.imshow(dst, clim=(vmin, vmax), cmap=cmap, interpolation="nearest", aspect='auto')
    plt.colorbar()

    ynum = dst.shape[0]/2
    if (isinstance(tads, pd.DataFrame)):
        for tad in tads.itertuples(name=None):
            x1 = tad[2]/resolution - s
            x2 = tad[3]/resolution - s
            x1 *= 1.41
            x2 *= 1.41
            if (x1 >0):
                xmed = (x1 + min([x2, (e-s)*1.41]))/2
                plt.plot([x1, min([xmed, x2])],[ynum, ynum-(xmed-x1)],
                         color='k', linestyle='dashed', linewidth=0.6)
            if (x2 < (e-s)*1.41):
                xmed = (max([x1, 0]) + x2) /2
                plt.plot([x2, xmed],[ynum, ynum-(x2-xmed)],
                         color='k', linestyle='dashed', linewidth=0.6)

    if (isinstance(loops, pd.DataFrame)):
        for loop in loops.itertuples(name=None):
            x1 = (loop[2] + loop[3])/2/resolution - s
            x2 = (loop[5] + loop[6])/2/resolution - s
            x1 *= 1.41
            x2 *= 1.41
            if (x1 >0):
                xmed = (x1 + min([x2, (e-s)*1.41]))/2
                plt.plot([x1, min([xmed, x2])],[ynum, ynum-(xmed-x1)],
                         color='b', linestyle='dashed', linewidth=0.5)
            if (x2 < (e-s)*1.41):
                xmed = (max([x1, 0]) + x2) /2
                plt.plot([x2, xmed],[ynum, ynum-(x2-xmed)],
                         color='b', linestyle='dashed', linewidth=0.6)

    if (distancemax > 0):
        ystart = max(0, ynum - distancemax/resolution)
        plt.ylim(ynum, ystart)
        pltyticks(ynum, ystart, 0, min(distancemax, (figend - figstart)), 5)
    else:
        ystart = 0
        plt.ylim(ynum, ystart)
        pltyticks(ynum, ystart, 0, figend - figstart, 5)

    if (label != ""): plt.title(label)
    if (xticks):
        pltxticks(0, (e-s)*1.41, figstart, figend, 10)
    else:
        xtickoff(plt)
