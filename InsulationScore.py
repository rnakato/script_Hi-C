#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from loadData import loadDenseMatrix

def calceach(mat, squaresize, resolution):
    matsize = int(squaresize / resolution)
    array = np.zeros(mat.shape[0])

    for i in range(mat.shape[0]):
        if(i - matsize < 0 or i + matsize >= mat.shape[0]): continue
        array[i] = mat[i-matsize:i, i+1:i+matsize+1].mean()
#    print(squaresize, np.nanmean(array))

    if (np.nanmean(array) != 0):
        array = np.log1p(array/np.nanmean(array))

    return array

def calcInsulationScore(mat, max_sqsize, step, resolution):
    imax = int(max_sqsize/step)
    for i in range(imax, 0, -1):
        if i==imax:
            InsulationScore = calceach(mat, i * step, resolution)
        else:
            InsulationScore = np.c_[InsulationScore, calceach(mat, i * step, resolution)]

    InsulationScore = InsulationScore.T
    df = pd.DataFrame(InsulationScore)
    df.index = np.arange(imax, 0, -1) * step
    df.columns = df.columns * resolution
    return df

class MultiInsulationScore:
    def __init__(self, mat, max_sqsize, step, resolution):
        imax = int(max_sqsize/step)
        for i in range(imax, 0, -1):
            if i==imax:
                InsulationScore = calceach(mat, i * step, resolution)
            else:
                InsulationScore = np.c_[InsulationScore, calceach(mat, i * step, resolution)]

        InsulationScore = InsulationScore.T
        InsulationScore = np.nan_to_num(InsulationScore)  # fill 0 in NA

        self.MI = pd.DataFrame(InsulationScore)
        self.MI.index = np.arange(imax, 0, -1) * step
        self.MI.columns = self.MI.columns * resolution

    def getInsulationScore(self, *, distance=500000):
        i = np.where(self.MI.index == distance)[0][0]
        return self.MI.iloc[i:i+1].T
    def getMultiInsulationScore(self):
        return self.MI.T

def getTADboundary(array, resolution):
    distance = int(100000 / resolution)
    slop = np.zeros(array.shape[0])
    for i in range(distance, array.shape[0] - distance):
        slop[i] = array[i - distance] - array[i + distance]

    boundary = []
    for i in range(1, len(slop)):
        if(slop[i-1] > 0 and slop[i] < 0 and array[i] <= -0.1):
            boundary.append(i)
    return boundary

def getInsulationScoreOfMultiSample(samples, labels, *, distance=500000):
    for i, sample in enumerate(samples):
        if i==0: Matrix = sample.getInsulationScore(distance=distance)
        else:    Matrix = pd.concat([Matrix, sample.getInsulationScore(distance=distance)], axis=1)
    Matrix = Matrix.replace(-np.inf,np.nan).fillna(0)
    Matrix.columns = labels
    return Matrix

def getDiffInsulationScoreOfMultiSample(samples, labels, *, distance=500000):
    for i, sample in enumerate(samples):
        df = sample.getInsulationScore(distance=distance) - samples[0].getInsulationScore(distance=distance)
        if i==0: pass
        elif i==1: DMat = df
        else:      DMat = pd.concat([DMat, df], axis=1)
    DMat = DMat.replace(-np.inf,np.nan).fillna(0)
    DMat.columns = labels
    return DMat

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", help="Input matrix", type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="Chromosome", type=str)
    parser.add_argument("resolution", help="Resolution of the input matrix", type=int)
#    parser.add_argument("--maxdistance", help="Max distance of MultiInsulationScore (default: 1000000)", type=int, default=1000000)
#    parser.add_argument("--step", help="stepsize of MultiInsulationScore (default: 100000)", type=int, default=100000)
    parser.add_argument("--num4norm", help="Read number after normalization (default: 10000000)", type=int, default=10000000)
    parser.add_argument("--distance", help="Distance of Insulation Score (default: 500000)", type=int, default=500000)

    args = parser.parse_args()
#    print(args)

#    if args.step < args.resolution:
#        print("Error: please set step (" + str(args.step) + ") larger than resolution (" + str(args.resolution) + ")")
#        sys.exit()

#    if args.distance > args.maxdistance:
#        print("Error: please set distance (" + str(args.distance) + ") less than maxdistance (" + str(args.maxdistance) + ")")
#        sys.exit()

    matrix = loadDenseMatrix(args.matrix)
    matrix = matrix * args.num4norm / np.nansum(matrix)

#    MI = MultiInsulationScore(matrix.values, args.maxdistance, args.step, args.resolution)
    MI = MultiInsulationScore(matrix.values, args.distance*2, args.distance, args.resolution)

#    print(MI.getMultiInsulationScore())

    # output InsulationScore to BedGraph
    df = MI.getInsulationScore(distance=args.distance)

    df = df.replace([np.inf, -np.inf], 0)

#    pd.set_option('display.max_rows', 1000)
#    print(df)
#    sys.exit()

    df.columns = ["Insulation Score"]
    df["chr"] = args.chr
    df["start"] = df.index
    df["end"] = df["start"] + args.resolution
    df = df.loc[:,["chr","start","end","Insulation Score"]]
    df.to_csv(args.output + ".bedGraph", sep="\t", header=False, index=False)
