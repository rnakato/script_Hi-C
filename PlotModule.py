import numpy as np
import matplotlib.pyplot as plt

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
