#!/usr/bin/env python
version="REPLACE_WITH_COMMIT_OR_VERSION"

import lib10x
import sys
import argparse
import numpy
import scipy
import time
import pybedtools
import pandas 
import os
import scipy.sparse
import scipy.stats
import gc

start_time=time.time()

parser = argparse.ArgumentParser()
parser.add_argument("gridFile1",help="1st set of grids BED formatted ")
parser.add_argument("fragFile",help="10x fragments BED formatted")
parser.add_argument("hdf5File",help="10x statistics HDF5 formatted")
parser.add_argument("outputFile",help="output prefix")
parser.add_argument("-b","--gridFile2",default="none",help="2nd set of grids BED formatted (optional)")
parser.add_argument("-e","--expectFalseDiscover",default=1000,help="expected false discoveries, [default: %default]")
parser.add_argument("-m","--maxTagFraction",default="none",help="max tag fraction to be considerrred, [default: %default]")
parser.add_argument("-g","--genomeSize",default=3000000000,help="bp size of the whole genome, [default: %default]")
parser.add_argument("-r","--regionSize",default=100000,help="bp size of the whole genome, [default: %default]")

try:
  args = parser.parse_args()
  gridFile1 = args.gridFile1
  gridFile2 = args.gridFile2
  fragFile = args.fragFile
  hdf5File = args.hdf5File
  outputFile = args.outputFile
  maxTagFraction = args.maxTagFraction
  expectFalseDiscover = args.expectFalseDiscover
  genomeSize = args.genomeSize
  regionSize = args.regionSize
except SystemExit:
  args = dict()
  gridFile1 = os.environ['HOME']+"/ws/hg/hg19/human_g1k_v37.ucsc.fasta.bed.w100Ks50K"
  gridFile2 = "none"
  maxTagFraction = "none"
  fragFile = "P04979_CR1.wg10.bam.fragBar.bed.a50000"
  hdf5File = "P04979_CR1.wg10.bam.a50000.h5"
  outputFile = "P04979_CR1.wg10.bam.grid100K.txt"
  genomeSize = 3000000000   #size of genome 3G
  regionSize = 100000       #size of the region 100K
  expectFalseDiscover = 1000 #tolerated false discoveries

print >>sys.stderr, "-Info: vesion:", version
print >>sys.stderr, "-Info: invoking command:", args
print >>sys.stderr, "-Info: gridFile1:", gridFile1
print >>sys.stderr, "-Info: gridFile2:", gridFile2
print >>sys.stderr, "-Info: fragFile:", fragFile
print >>sys.stderr, "-Info: hdf5File:", hdf5File
print >>sys.stderr, "-Info: outputFile:", outputFile
print >>sys.stderr, "-Info: genomeSize:", genomeSize
print >>sys.stderr, "-Info: regionSize:", regionSize
print >>sys.stderr, "-Info: expectFalseDiscover:", expectFalseDiscover

#NOTE:
#  big memory objects: 
#    overlap1, overlap2, fragdata, frag2tag
#    X_i, X_j, X_data, X
#    Y_j, Y_k, Y_data, Y
gridLabels = ["chr","start","end"]
grid1 = pybedtools.BedTool(gridFile1) #.to_dataframe(comment='#',names=gridLabels,low_memory=False,sep="\t")
print >>sys.stderr, "grid1Size:", len(grid1)
overlapFile1 = ".".join([outputFile,"overlap1","mtx","gz"])
overlapCMD1 = "gridOverlap.R %s %s %s" % (gridFile1, fragFile, overlapFile1)
print >>sys.stderr, overlapCMD1
os.system(overlapCMD1)
overlap1 = pandas.read_csv(overlapFile1,delimiter='\t')
print >>sys.stderr, "overlap1Size:", len(overlap1)
if gridFile2 == "none":
  #overlap2 = overlap1
  overlap2 = None
  print >>sys.stderr, "overlap2 = overlap1"
else:
  grid2 = pybedtools.BedTool(gridFile2) #.to_dataframe(comment='#',names=gridLabels,low_memory=False,sep="\t")
  print >>sys.stderr, "grid2Size:", len(grid2) 
  overlapFile2 = ".".join([outputFile,"overlap2","mtx","gz"])
  overlapCMD2 = "gridOverlap.R %s %s %s" % (gridFile2, fragFile, overlapFile2)
  print >>sys.stderr, overlapCMD2
  os.system(overlapCMD2)
  overlap2 = pandas.read_csv(overlapFile2,delimiter='\t')
  print >>sys.stderr, "overlap2Size:", len(overlap2)

fragLabels = ["chr","start","end","fragTagID","fragSize","tagFragFreq"]
#fragdata = pybedtools.BedTool(fragFile).to_dataframe(comment='#',names=fragLabels,low_memory=False,sep="\t")
fragdata = lib10x.bed2df(fragFile,fragLabels)
#assert len(fragdata) == len(pybedtools.BedTool(fragFile).to_dataframe(comment='#',names=fragLabels,low_memory=False,sep="\t"))
print >>sys.stderr, "fragdataSize:", len(fragdata)

#ASSUME:
#  hdf5 to contain ALL possible Tb's
#  fragFile to caontain ALL Fragments
regionFraction = regionSize/float(genomeSize)              #fraction of the size of a region compared to a Gennome
print >>sys.stderr, "regionFraction:", regionFraction
tagFractions = pandas.read_hdf(hdf5File,key="T")
print >>sys.stderr, "tagFractions:", tagFractions
fragNumber = pandas.read_hdf(hdf5File,key="N_f").iloc[0,0]
print >>sys.stderr, "fragNumber:", fragNumber
averageFragLength = pandas.read_hdf(hdf5File,key="L_f").iloc[0,0] #size of a fragment in bases
print >>sys.stderr, "averageFragLength:", averageFragLength
tagFractionSquareSum=numpy.sum(tagFractions*tagFractions,axis=1)[0]
print >>sys.stderr, "tagFractionSquareSum:", tagFractionSquareSum
expectShareFrag=(fragNumber**2)*(regionFraction**2)*tagFractionSquareSum #FIX, u0, expected shared Fragment by random chance
print >>sys.stderr, "expectShareFrag:", expectShareFrag
expectJunctionFrag=fragNumber*averageFragLength/genomeSize                    #u, expected shared Fragment by junction formation
print >>sys.stderr, "expectJunctionFrag:", expectJunctionFrag

if maxTagFraction == "none":
  maxTagFraction = max(numpy.percentile(tagFractions,99),0.001) 
else:
  maxTagFraction = float(maxTagFraction)

#OLD
#print >>sys.stderr, "creating hash table for tag and fraction..."
#gtk1_time = time.time()
#tidx2tag = [-1] * len(tagFractions.iloc[0]) # 0-based, from filtered idx to original idx
#tag2tidx = [-1] * len(tagFractions.iloc[0]) # 0-based, from original idx to filtered idx
#tidx = 0
#for i in range(0, len(tidx2tag)):
#  if tagFractions.iloc[0,i] <= maxTagFraction:
#    tidx2tag[tidx] = i
#    tag2tidx[i] = tidx
#    tidx = tidx + 1
#
#tidx2tag = tidx2tag[:tidx]
#frag2tag = list(fragdata.iloc[:,3]-1)
#gtk2_time = time.time()
#print >>sys.stderr, gtk2_time-gtk1_time, "s done"

#IMPROVE: 
#  Get a set of overaboundant tags from tagFractions
#  Replace X.iloc[:,t] by 0
excludeTag = set([ ti for ti in xrange(0, len(tagFractions.iloc[0,:])) if tagFractions.iloc[0,ti] < maxTagFraction ])

print >>sys.stderr, "creating matrix X..."
#DONE: X(i,j) is overlaps of Region1(i) and Fragment with Tb(j)
gtk1_time = time.time()
#X_i = list(overlap1.iloc[:,0]-1)   #convert to 0-based python index
#X_j = [tag2tidx[frag2tag[i]] for i in overlap1.iloc[:,1]-1]
#X_j = list(overlap1.iloc[:,1]-1)
#X_data = X_data[ X_data['j'] >= 0 ]
X_i = [ overlap1.iloc[oi,0] for oi in xrange(0, len(overlap1)) if overlap1.iloc[oi,1] not in excludeTag ]
X_j = [ fragdata['fragTagID'][overlap1.iloc[oi,1]-1]-1 for oi in xrange(0, len(overlap1)) if overlap1.iloc[oi,1] not in excludeTag ]
X_data = pandas.DataFrame({ 'i' : X_i, 'j' : X_j }) # region-i, tag-j
X_data = X_data.groupby(['i', 'j']).size()
X_data = X_data.reset_index()
X_data.columns = ['i', 'j', 'v']
X_shape = (len(grid1), len(tagFractions.iloc[0,:]))
X = scipy.sparse.coo_matrix((X_data['v'], (X_data['i'], X_data['j'])), shape=X_shape)
X = X.tocsr()
del X_i, X_j, X_data
gc.collect()
print >>sys.stderr, X.shape
print >>sys.stderr, X.nnz
gtk2_time = time.time()
print >>sys.stderr, gtk2_time-gtk1_time, "s done"

print >>sys.stderr, "creating matrix Y..."
#DONE: Y(j,k) is indicator if overlap of Region2(k) and Fragment with Tb(j)
gtk1_time = time.time()
if gridFile2 == "none":
  Y = X.transpose(copy=True)
else:
  Y_k = [ overlap2.iloc[oi,0] for oi in xrange(0, len(overlap2)) if overlap2.iloc[oi,1] not in excludeTag ]
  Y_j = [ fragdata['fragTagID'][overlap2.iloc[oi,1]-1]-1 for oi in xrange(0, len(overlap1)) if overlap2.iloc[oi,1] not in excludeTag ]
  #Y_k = [ overlap2.iloc[oi,1] for oi in xrange(0, len(overlap2)) if overlap1.iloc[oi,1] not in excludeTag ]
  #Y_j = [tag2tidx[frag2tag[i]] for i in overlap1.iloc[:,1]-1]
  #Y_k = list(overlap1.iloc[:,0]-1)   #convert to 0-based python index
  Y_data = pandas.DataFrame({ 'j' : Y_j, 'k' : Y_k }) # region-i, tag-j
  #Y_data = Y_data[ Y_data['j'] >= 0 ]
  Y_data = Y_data.groupby(['j', 'k']).size()
  Y_data = Y_data.reset_index()
  Y_data.columns = ['j', 'k', 'v']
  Y_shape = (len(tidx2tag), len(grid2))
  Y = scipy.sparse.coo_matrix((Y_data['v'], (Y_data['j'], X_data['k'])), shape=Y_shape)
  Y = Y.tocsc()
  del Y_j, Y_k, Y_data
  gc.collect()

Y = Y>0
print >>sys.stderr, Y.shape
print >>sys.stderr, Y.nnz
gtk2_time = time.time()
print >>sys.stderr, gtk2_time-gtk1_time, "s done"

#DONE create a sparse matrix M = X * t(Y)
print >>sys.stderr, "multiplying matrix X and Y..."
gtk1_time = time.time()
M = X.dot(Y)
del X, Y, overlap1, overlap2, fragdata
gc.collect()
print >>sys.stderr, M.shape
print >>sys.stderr, M.nnz
gtk2_time = time.time()
print >>sys.stderr, gtk2_time-gtk1_time, "s done"

#DEBUG: 
#  test computation limitations
#  X = scipy.sparse.rand(60000, 60000, density=0.002, format='csr')
#  Y = X.transpose(copy=True)
#  M = X.dot(Y)   #memory bust

#compoute tcut from FDR control: affordable computations N_c = M^2 * (1-ppois(u0,t))

cutFDR = scipy.stats.poisson.ppf(1-expectFalseDiscover/float(M.shape[0]**2),expectShareFrag)
print >>sys.stderr, "cutFDR:", cutFDR
M = M>cutFDR
print >>sys.stderr, M.nnz
M = M.tocoo()
M_data = pandas.DataFrame({ 'grid1' : M.row, 'grid2' : M.col, 'm12' : M.data })
M_data = M_data[ M_data['m12'] > cutFDR ]
M_data = M_data.sort_values(by=['m12','grid1','grid2'], ascending=[0,1,1])
print >>sys.stderr, "writing results..."
M_data.to_csv(".".join([outputFile,"RData"]), sep="\t")
end_time=time.time()
print >>sys.stderr, "total", end_time-start_time, "s"
print >>sys.stderr, "GRIDSCANDONE"
