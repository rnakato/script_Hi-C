#! /usr/bin/env python
# -*- coding: utf-8 -*- 

import numpy as np
import pandas as pd
import argparse
from loadData import loadJuicerMatrix

def calcDI(mat, resolution):
    def getDI(mat, i, len):
        A = np.triu(mat[i-len:i, i-len:i]).sum()
        B = np.triu(mat[i:i+len, i:i+len]).sum()
        E = (A + B)/2
        temp = np.nan_to_num(((A-E)**2)/E) + np.nan_to_num(((B-E)**2)/E)
        DI = np.nan_to_num((B-A)/np.abs(B-A)) * temp
        return DI
    
    len = int(2000000 / resolution)
    array = np.zeros(mat.shape[0])
    for i in range(len, mat.shape[0]-len):
        array[i] = getDI(mat, i, len)
    return array

def getDifferentialIndexOfMultiSample(samples, labels):
    for i, sample in enumerate(samples):
        if i==0: Matrix = sample.DI
        else:    Matrix = np.vstack((Matrix,sample.DI)) 
    Matrix = pd.DataFrame(Matrix)
    Matrix = Matrix.replace(-np.inf,np.nan).fillna(0)
    Matrix.index = labels
    return Matrix

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", help="Input matrix", type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="Chromosome", type=str)
    parser.add_argument("resolution", help="Resolution of the input matrix", type=int)
    parser.add_argument("--num4norm", help="Read number after normalization (default: 10000000)", type=int, default=10000000)
    parser.add_argument("--distance", help="Distance of Insulation Score (default: 500000)", type=int, default=500000)

    args = parser.parse_args()
    print(args)

    matrix = loadJuicerMatrix(args.matrix)
    matrix = matrix * args.num4norm / np.nansum(matrix)

    MI = MultiInsulationScore(matrix.values, 1000000, 100000, args.resolution)

    # output InsulationScore to BedGraph
    df = MI.getInsulationScore(distance=args.distance)
    df = df.replace([np.inf, -np.inf], 0)
    df.columns = ["Insulation Score"]
    df["chr"] = args.chr
    df["start"] = df.index
    df["end"] = df["start"] + args.resolution
    df = df.loc[:,["chr","start","end","Insulation Score"]]
    df.to_csv(args.output + ".bedGraph", sep="\t", header=False, index=False)
