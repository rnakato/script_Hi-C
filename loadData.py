import pandas as pd

def loadDenseMatrix(filename):
    print(filename)
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    return data

def loadTADs(filename, chr, *, start=0, end=99999999999):
    tads = pd.read_csv(filename, delimiter='\t', usecols=['chr1','x1','x2'])
    tads = tads[tads.chr1 == chr]
    tads = tads[tads.x1 < end]
    tads = tads[tads.x2 >= start]
    return tads

def readBedGraph(file, chr, *, start=-1, end=-1):
    bedgraph = pd.read_csv(file, delimiter='\s+', header=None)
    bedgraph.columns = ["chr", "start", "end", "value"]
    bedgraph = bedgraph[bedgraph["chr"] == chr]
    bedgraph = bedgraph.set_index("start")
    bedgraph = bedgraph["value"]

    if start >= 0 and end > 0:
        bedgraph = bedgraph[bedgraph.index >= start]
        bedgraph = bedgraph[bedgraph.index <= end]
    return bedgraph
