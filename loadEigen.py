def loadEigen(filename, refFlat, chr, res):
   print(filename)
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
このように関数を定義して、
eigen = loadEigen("/media/DATAPART3/Hi-C/Sakata/RPE/Juicer/2017_009A_R_Ct/matrix/intrachromosomal/50000/eigen.KR.chr21.txt.gz",
        "/home/Database/UCSC/hg19/refFlat.dupremoved.txt",
        "chr21",
        50000)
