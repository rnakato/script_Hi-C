import numpy as np
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

def axxticks(ax, start, end, figstart, figend, nmem):
    mem = int((end - start)/nmem)
    x = start + np.arange(nmem+1) * mem
    val = (figend - figstart) / nmem
    xval = (figstart + np.arange(nmem+1)*val)/Mega
    ax.set_xticks(x)
    ax.set_xticklabels(xval)

def xtickoff(plt):
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False  # labels along the bottom edge are off
    )

def drawHeatmapSquare(plt, matrix, resolution, *, tads="", figstart=0, figend=0,
                      vmin=0, vmax=50, label="", xticks=True):
    s = int(figstart / resolution)
    e = int(figend   / resolution)
    if (e==0): e = matrix.shape[0]
    plt.imshow(matrix.iloc[s:e,s:e], clim=(vmin, vmax),
               cmap=generate_cmap(['#FFFFFF', '#d10a3f']),
               interpolation="nearest", aspect='auto')

    if (type(tads) != "str"):
        for tad in tads.itertuples(name=None):
            x1 = int(int(tad[2])/resolution) - s
            x2 = int(int(tad[3])/resolution) - s
            plt.vlines(x1, 0, matrix.shape[0], linestyle='dashed', linewidth=0.3)
            plt.vlines(x2, 0, matrix.shape[0], linestyle='dashed', linewidth=0.3)

    if (label != ""): plt.title(label)
    if (xticks):
        pltxticks(0, e-s, figstart, figend, 10)
    else:
        xtickoff(plt)

def drawHeatmapTriangle(plt, matrix, resolution, *, tads="", figstart=0, figend=0,
                        vmin=0, vmax=50, label="", xticks=True):
    s = int(figstart / resolution)
    e = int(figend   / resolution)
    if (e==0): e = matrix.shape[0]
    dst = ndimage.rotate(matrix.iloc[s:e,s:e], 45,
                         order=0, reshape=True, prefilter=False, cval=0)
    plt.imshow(dst, clim=(vmin, vmax), cmap=generate_cmap(['#FFFFFF', '#d10a3f']),
               interpolation="nearest", aspect='auto')

    ynum = dst.shape[0]/2
    if (type(tads) != "str"):
        for tad in tads.itertuples(name=None):
            x1 = (int(tad[2])/resolution) - s
            x2 = (int(tad[3])/resolution) - s
            x1 *= 1.41
            x2 *= 1.41
            xmed = (x1+x2)/2
            plt.plot([x1, xmed],[ynum, ynum-(xmed-x1)], color='k', linestyle='dashed', linewidth=0.6)
            plt.plot([x2, xmed],[ynum, ynum-(xmed-x1)], color='k', linestyle='dashed', linewidth=0.6)

    plt.ylim(ynum,0)
    if (label != ""): plt.title(label)
    if (xticks):
        pltxticks(0, (e-s)*1.41, figstart, figend, 10)
    else:
        xtickoff(plt)
