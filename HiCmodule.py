import numpy as np
import scipy.stats as sp
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from InsulationScore import *
from DirectionalityIndex import *
from loadData import loadJuicerMatrix
from generateCmap import generate_cmap

def getNonZeroMatrix(A, lim_pzero):
    A = A.fillna(0)
    pzero_row = A[A>0].count(axis=0)/A.shape[0]
    index = pzero_row[pzero_row > lim_pzero].index
    pzero_col = A[A>0].count(axis=1)/A.shape[1]
    columns = pzero_col[pzero_col > lim_pzero].index

    A = A[index]
    A = A.loc[columns]

    return A

def loadEigen(filename, refFlat, chr, res):
    print(filename)
    if filename == "": return
    
    eigen = np.loadtxt(filename)
    gene = pd.read_csv(refFlat, delimiter='\t', header=None, index_col=0)
    gene = gene[gene.iloc[:,1] == chr]
    gene = gene.iloc[:,5].values

    genenum = np.zeros(len(eigen), int)
    for row in gene:
        if int(row/res)>= len(eigen): continue
        genenum[int(row/res)] += 1

    from scipy.stats import pearsonr
    if pearsonr(genenum[~np.isnan(eigen)], eigen[~np.isnan(eigen)])[0] < 0:
        eigen = -eigen
        
    return eigen

class JuicerMatrix:
    def __init__(self, norm, rawmatrix, oematrix, eigenfile, refFlat, chr, res):
        self.res = res
        self.raw = loadJuicerMatrix(rawmatrix)
        self.oe  = loadJuicerMatrix(oematrix)
        self.eigen = loadEigen(eigenfile, refFlat, chr, res)
        if norm == "RPM":
            self.raw = self.raw * 10000000 / np.nansum(self.raw)
            self.oe  = self.oe  * 10000000 / np.nansum(self.oe)
        self.InsulationScore = MultiInsulationScore(self.getmatrix().values, 1000000, 100000, self.res)

    def getmatrix(self, *, isOE=False, isNonZero=False):
        if isOE == False:
            if isNonZero == True:
                return getNonZeroMatrix(self.raw, 0)
            else:
                return self.raw
        else:
            if isNonZero == True:
                return getNonZeroMatrix(self.oe, 0)
            else:
                return self.oe
        
    def getlog(self, *, isOE=False, isNonZero=False):
        mat = self.getmatrix(isOE=isOE, isNonZero=isNonZero)
        logmat = mat.apply(np.log1p)
        return logmat

    def getZscore(self, *, isOE=False):
        logmat = self.getlog(isOE=isOE)
        zmat = pd.DataFrame(sp.stats.zscore(logmat, axis=1), index=logmat.index, columns=logmat.index)
        return zmat

    def getPearson(self, *, isOE=False):
        oe = self.getlog(isOE=isOE)
#        oe = self.getmatrix(isOE=isOE)
#        cov = logmat.cov()
 #       ccmat = np.corrcoef(cov)
        ccmat = np.corrcoef(oe)
        ccmat[np.isnan(ccmat)] = 0
        ccmat = pd.DataFrame(ccmat, index=oe.index, columns=oe.index)
        return ccmat
    
    def getEigen(self):
        from sklearn.decomposition import PCA
        pca = PCA()
        ccmat = self.getPearson()
#        index = np.isnan(ccmat).all(axis=1)
        transformed = pca.fit_transform(ccmat)
        pc1 = transformed[:, 0]
#        pc1[index] = np.nan
        return transformed[:, 0]

    def getInsulationScore(self, *, distance=500000):
        return self.InsulationScore.getInsulationScore(distance=distance)

    def getMultiInsulationScore(self):
        return self.InsulationScore

    def getDirectionalityIndex(self, *, distance=1000000):
        return calcDI(self.raw.values, self.res, distance=distance)

def ExtractMatrix(mat,s,e):
    if e==-1:
        return mat.values[s:,s:]
    else:
        return mat.values[s:e,s:e]
    
def ExtractMatrixIndex(mat,index1,index2):
    mat = mat[index1,:]
    mat = mat[:,index2]
    return mat

def ExtractTopOfPC1(mat, pc1, nbin):
    sortedindex = np.argsort(pc1)
    sortmat = mat[sortedindex,:]
    sortmat = sortmat[:,sortedindex]
    sortmat = np.concatenate((sortmat[:nbin,:] , sortmat[-nbin:,:]), axis=0)
    sortmat = np.concatenate((sortmat[:,:nbin] , sortmat[:,-nbin:]), axis=1)
    return sortmat

def getCompartment(mat, pc1):
    indexA = pc1 > 1
    indexB = pc1 < -1
    A = mat[indexA]
    A = A[:,indexA]
    B = mat[indexB]
    B = B[:,indexB]
    return A, B

def getCompartment_inter(mat, pc1_odd, pc1_even, nbin):
    sortedindex_odd = np.argsort(pc1_odd)
    sortedindex_even = np.argsort(pc1_even)
    sortmat = mat[sortedindex_odd,:]
    sortmat = sortmat[:,sortedindex_even]
    A = sortmat[:nbin,:nbin]
    B = sortmat[-nbin:,-nbin:]
#    indexA_odd = pc1_odd > 1
#    indexB_odd = pc1_odd < -1
 #   indexA_even = pc1_even > 1
  #  indexB_even = pc1_even < -1
  #  A = mat[indexA_odd]
  #  A = A[:,indexA_even]
  #  B = mat[indexB_odd]
   # B = B[:,indexB_even]
    return A, B
