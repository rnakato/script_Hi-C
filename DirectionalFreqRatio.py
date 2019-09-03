import numpy as np

def make3dmatrixRatio(samples, smoooth):
    from scipy import ndimage
    n = len(samples)
    Ct = ndimage.median_filter(samples[0].getlog(isNonZero=False), smoooth)
    x, y = Ct.shape
    for i, sample in enumerate(samples[1:]):
        if i==0:
            data = sample.getlog(isNonZero=False)
            Matrix = ndimage.median_filter(data - Ct, smoooth)
        else:
            data = sample.getlog(isNonZero=False)
            M = ndimage.median_filter(data - Ct, smoooth)
            Matrix = np.concatenate((Matrix, M))
    Matrix = Matrix.reshape(n-1,x,y)
    return Matrix

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

class DirectionalFreqRatio:
    def __init__(self, mat, resolution, *, distance=2000000):
        self.arrayplus = getDirectionalFreqRatio(mat, resolution, "+", distance=distance)
        self.arrayminus = getDirectionalFreqRatio(mat, resolution, "-", distance=distance)

    def getarrayplus(self):
        return self.arrayplus

    def getarrayminus(self):
        return self.arrayminus

    def getarraydiff(self):
        return self.arrayplus - self.arrayminus
