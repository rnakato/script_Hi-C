import pandas as pd

def readBedGraph(file, chr):   
    bedgraph = pd.read_csv(file, skiprows=5, delimiter=' ', header=None)
    bedgraph.columns = ["chr", "start", "end", "value"]
    bedgraph = bedgraph[bedgraph["chr"] == chr]
    bedgraph = bedgraph.set_index("start")
    bedgraph = bedgraph["value"]
    return bedgraph

