#!/usr/bin/env python
#using two spaces as indent

import sys
import argparse
import numpy
import scipy
import time
import pandas 
import os
import scipy.sparse
import scipy.stats
import scipy.io
import gc
import shelve
import tempfile
import json

try:
  import lib10x
except ImportError:
  import zoomx.lib10x
start_time=time.time()

def main(args):

  """
    GridScan Module of ZoomX for Poisson Scan of 10X Genomics Bam

    INPUT:
      gridFile1 (1st set of grids BED formatted):  
        $hg19.bed.w10Ks5K
      fragFile (10x fragments BED formatted):
        $bam.multiFragBar.bed.a0
      hdf5File (10x statistics HDF5 formatted):
        $bam.multiFragBar.h5.a0
      outputPrefix (output prefix):
        $bam
        

    VARIABLE:

    OUTPUT:

  """

  if args:
    gridFile1 = args.gridFile1
    gridFile2 = args.gridFile2
    skipFile1 = args.skipFile1
    skipFile2 = args.skipFile2
    fragFile = args.fragFile
    hdf5File = args.hdf5File
    outputPrefix = args.outputPrefix
    maxTagFraction = args.maxTagFraction
    expectFalseDiscover = args.expectFalseDiscover
    aneuploidyAllowance = args.aneuploidyAllowance
    maxDuplication = args.maxDuplication
    maxJunction = args.maxJunction
    genomeSize = args.genomeSize
    gridSize = args.gridSize
    numSub = args.numSub
    tempdir = args.tempdir
    debug = args.debug
  else:
    bamfile="LS411N_CR1.wg12.bam"
    if not os.path.isfile(bamfile): quit() # quit if not in test directory
    print >>sys.stderr, "fall back to run program tests..."
    gridFile1 = "human_g1k_v37.fasta.bed.w10Ks5K"
    gridFile2 = "same"
    skipFile1 = "none"
    #skipFile1 = os.environ['HOME']+"/ws/hg/hg19/human_g1k_v37.fasta.acen+gap+hc.bed"
    skipFile2 = "same"
    maxTagFraction = 0.99
    fragFile = "LS411N_CR1.wg12.bam.multiFragBar.bed.a0"
    hdf5File = "LS411N_CR1.wg12.bam.multiFragBar.h5.a0"
    outputPrefix = bamfile 
    genomeSize = 1000000   #size of genome region of bam
    numSub = 3
    gridSize = 10000          #size of the region 10K
    expectFalseDiscover = 10  #tolerated false discoveries
    maxDuplication = 1000
    maxJunction = 1000
    aneuploidyAllowance = 4
    tempdir = 'zoomx_'+lib10x.randascii(4)
    debug = False
  
  lib10x.setTemporaryDir(tempdir)
  scanConfig = dict()
  
  #print >>sys.stderr, "-Info: invoking arguments:", args
  sys.stderr.write("gridScan.py {} {} {} {}\n".format(gridFile1,fragFile,hdf5File,outputPrefix))
  sys.stderr.write("args="+str(args)+"\n")
  if debug:
    print >>sys.stderr, "-Info: gridFile1:", gridFile1
    print >>sys.stderr, "-Info: gridFile2:", gridFile2
    print >>sys.stderr, "-Info: skipFile1:", skipFile1
    print >>sys.stderr, "-Info: skipFile2:", skipFile2
    print >>sys.stderr, "-Info: fragFile:", fragFile
    print >>sys.stderr, "-Info: hdf5File:", hdf5File
    print >>sys.stderr, "-Info: outputPrefix:", outputPrefix
    print >>sys.stderr, "-Info: genomeSize:", genomeSize
    print >>sys.stderr, "-Info: gridSize:", gridSize
    print >>sys.stderr, "-Info: expectFalseDiscover:", expectFalseDiscover
  
  #NOTE:
  #  big memory objects: 
  #    overlap1, overlap2, fragdata, frag2tagidx
  #    X_i, X_j, X_data, X
  #    Y_j, Y_k, Y_data, Y
  overlapLabels = ["chr","start","end","overlaps"]
  overlapDtypes = dict(zip(overlapLabels, [ object, pandas.np.int64, pandas.np.int64, pandas.np.int64 ]))
  gridLabels = ["chr","start","end"]
  gridDtypes = dict(zip(gridLabels,["string",numpy.int32,numpy.int32]))
  bedpeLabels = ["chr1","start1","end1","chr2","start2","end2","m12"]
  fragLabels = ["chr","start","end","moltagidx","tagmolcnt","MI","molsize","mollen","HP","PS","BX"]
  fragDtypes = {"chr":'object',"start":'object',"end":'object',"moltagidx":numpy.int32,"tagmolcnt":numpy.int32,"MI":numpy.int32,"molsize":numpy.int32,"mollen":numpy.int32,"HP":numpy.int32,"PS":numpy.int32,"BX":'object'}
  
  frag = pandas.read_table(fragFile, names=fragLabels, dtype=fragDtypes)
  frag['index'] = xrange(1,frag.shape[0]+1)
  tmpFragFile = lib10x.getTemporaryFile('frag_bed') #1-based 
  frag.iloc[:,range(0,3)+[frag.shape[1]-1]].to_csv(tmpFragFile, header=False, sep="\t", index=False)
  if debug: 
    sys.stderr.write("tmpFragFile='"+tmpFragFile+"'\n")
  
  grid1 = pandas.read_table(gridFile1, names=gridLabels, dtype=gridDtypes)
  grid1['index'] = xrange(1,grid1.shape[0]+1)
  tmpGridFile1 = lib10x.getTemporaryFile('grid1_bed') #1-based 
  if debug:
    sys.stderr.write("tmpGridFile='"+tmpGridFile1+"'\n")
  grid1.to_csv(tmpGridFile1,header=False, sep="\t", index=False)
  grid1 = grid1.reindex()
  print >>sys.stderr, "grid1Size:", len(grid1)
    
  gtk1_time = time.time()
  
  overlapFile1 = lib10x.getTemporaryFile('overlap1_mtx_gz')
  if debug: 
    sys.stderr.write("overlapFile1='"+overlapFile1+"'\n") #1-based beda to bedb intersect file
  overlapCMD1 = "bedtools intersect -wa -wb -a %s -b %s | cut -f4,8 >%s" % (tmpGridFile1, tmpFragFile, overlapFile1)
  if debug: 
    sys.stderr.write("overlapCMD1="+str(overlapCMD1)+"\n")
    gtk1_time = time.time()
  if lib10x.run_cmd(overlapCMD1) == False:
    sys.stderr.write("Error in running overlapCMD1 overlapping..."); quit()
  if debug:
    gtk2_time = time.time()
    print >>sys.stderr, gtk2_time-gtk1_time, "s done"
  
  overlap1 = pandas.read_csv(overlapFile1,header=0,names=['grid','mol'],\
    dtype={'grid':pandas.np.int32,'mol':pandas.np.int32}, delimiter='\t')
  print >>sys.stderr, "overlap1Size:", len(overlap1)
  
  if gridFile2 == "same":
    #overlap2 = overlap1
    grid2 = grid1
    gridFile2 = gridFile1
    overlap2 = None
    print >>sys.stderr, "overlap2 = overlap1"
  else:
    grid2 = pandas.read_table(gridFile2, names=gridLabels, dtype=gridDtypes)
    grid2['index'] = xrange(0,grid2.shape[0])
    tmpGridFile2 = lib10x.getTemporaryFile('grid2_bed') #0-based 
    grid2.to_csv(tmpGridFile2,header=False, sep="\t", index=False)
    grid2 = grid2.reindex()
    print >>sys.stderr, "grid2Size:", len(grid2) 
    overlapFile2 = lib10x.getTemporaryFile('overlap2_mtx_gz') 
    overlapCMD2 = "bedtools intersect -wa -wb -a %s -b %s | cut -f4,8 >%s" % (tmpGridFile2, tmpFragFile, overlapFile2)
    if debug:
      sys.stderr.write("overlapCMD1="+str(overlapCMD2)+"\n")
      gtk1_time = time.time()
    if lib10x.run_cmd(overlapCMD2) == False:
      sys.stderr.write("Error in running overlapCMD2 overlapping..."); quit()
    if debug: 
      gtk2_time = time.time()
      print >>sys.stderr, gtk2_time-gtk1_time, "s done"
    overlap2 = pandas.read_csv(overlapFile2,header=0,names=['grid','mol'],\
  	  dtype={'grid':pandas.np.int32,'mol':pandas.np.int32}, delimiter='\t')
    print >>sys.stderr, "overlap2Size:", len(overlap2)
  
  #NOTE: make skipTable to be load
  #skipLabels = ["chr","start","end","count"]
  if skipFile1 != "none":
    skipOverlapFile1 = ".".join([outputPrefix,"skipOverlap1","txt","gz"])
    skipOverlapCMD1 = "bedtools intersect -c -a %s -b %s | gzip -c - > %s" % (gridFile1, skipFile1, skipOverlapFile1)
    if debug:
      print >>sys.stderr, "skipOverlapCMD1"
      print >>sys.stderr, skipOverlapCMD1
    if lib10x.run_cmd(skipOverlapCMD1) == False:
      sys.stderr.write("Error in running skipOverlapCMD1 overlapping..."); quit()
    skipOverlap1 = pandas.read_csv(skipOverlapFile1,header=None,names=overlapLabels,dtype=overlapDtypes,delimiter='\t')#,index_col=skipLabels)
    skipIndex1 = [ i for i in xrange(0,len(skipOverlap1)) if skipOverlap1.iloc[i,3]>0 ]
  else:
    skipIndex1 = []
  print >>sys.stderr, "skipIndex1 elements:", len(skipIndex1)
  
  if skipFile2 == "same":
    skipFile2 = skipFile1
    skipIndex2 = skipIndex1
  elif skipFile2 == "none":
    skipIndex2 = []
  else:
    skipOverlapFile2 = ".".join([outputPrefix,"skipOverlap2","txt","gz"])
    skipOverlapCMD2 = "bedtools intersect -c -a %s -b %s | gzip -c - > %s" % (gridFile2, skipFile2, skipOverlapFile2)
    if debug:
      print >>sys.stderr, "skipOverlapCMD2"
      print >>sys.stderr, skipOverlapCMD2
    os.system(skipOverlapCMD2)
    skipOverlap2 = pandas.read_csv(skipOverlapFile2,header=None,names=overlapLabels,dtype=overlapDtypes,delimiter='\t')#,index_col=skipLabels)
    skipIndex2 = [ i for i in xrange(0,len(skipOverlap2)) if skipOverlap2.iloc[i,3]>0 ]
  print >>sys.stderr, "skipIndex2 elements:", len(skipIndex2)
  
  #fragLabels = ["chr","start","end","fragTagID","fragSize","fragLength","tagFragFreq","MI","HP","PS"]
  #fragdata = pybedtools.BedTool(fragFile).to_dataframe(comment='#',names=fragLabels,low_memory=False,sep="\t")
  fragdata = lib10x.bed2df(fragFile,fragLabels)
  fragdataSize = len(fragdata)
  print >>sys.stderr, "fragdataSize:", fragdataSize
  #ASSUME:
  #  hdf5 to contain ALL possible Tb's
  #  fragFile to caontain ALL Fragments
  regionFraction = gridSize/float(genomeSize)              #fraction of the size of a region compared to a Gennome
  print >>sys.stderr, "regionFraction:", regionFraction
  tagFractions = pandas.read_hdf(hdf5File,key="T")
  print >>sys.stderr, "tagFractions:", tagFractions
  fragNumber = pandas.read_hdf(hdf5File,key="N_f").iloc[0,0]
  print >>sys.stderr, "fragNumber:", fragNumber
  averageFragLength = pandas.read_hdf(hdf5File,key="L_f").iloc[0,0] #size of a fragment in bases
  print >>sys.stderr, "averageFragLength:", averageFragLength
  tagFractionSquareSum=numpy.sum(tagFractions*tagFractions,axis=1)[0]
  print >>sys.stderr, "tagFractionSquareSum:", tagFractionSquareSum
  expectShareFrag=(fragNumber**2)*(regionFraction**2)*tagFractionSquareSum     #FIX, u0, expected shared Fragment by random chance
  print >>sys.stderr, "expectShareFrag:", expectShareFrag
  expectJunctionFrag=fragNumber*averageFragLength/float(genomeSize)                    #u, expected shared Fragment by junction formation
  print >>sys.stderr, "expectJunctionFrag:", expectJunctionFrag
  
  if maxTagFraction <0 or maxTagFraction>0.99:
    medTagFraction = numpy.median(tagFractions)
    madTagFraction = lib10x.mad(tagFractions)
    maxTagFraction = medTagFraction + 3*1.8*madTagFraction
    #maxTagFraction = max(numpy.percentile(tagFractions,99),0.0001) #don't cut off anything less than 1/10000 
  else:
    maxTagFraction = maxTagFraction
  
  maxTagFragments = maxTagFraction * fragNumber
  print >>sys.stderr, "maxTagFraction:", maxTagFraction
  print >>sys.stderr, "maxTagFragments:", maxTagFragments
  print >>sys.stderr, "per95TagFraction:", numpy.percentile(tagFractions,95)
  print >>sys.stderr, "per90TagFraction:", numpy.percentile(tagFractions,90)
  print >>sys.stderr, "creating hash table for tag and fraction..."

  tagcol2tagidx = [-1] * len(tagFractions.iloc[0]) # 0-based, from tag column idx in X to original tag idx in .h5
  tagidx2tagcol = [-1] * len(tagFractions.iloc[0]) # 0-based, from original tag idx in .h5 to tag column idx
  tcol = 0                               
  tagExcludeFile = ".".join([outputPrefix,"exclude","tagidx"])
  tagExcludeFileHandle = open(tagExcludeFile,"w")
  for i in range(0, len(tagcol2tagidx)):
    if tagFractions.iloc[0,i] <= maxTagFraction: #tag included
      tagcol2tagidx[tcol] = i  # col -> idx
      tagidx2tagcol[i] = tcol  # idx -> col 
      tcol = tcol + 1
    else:
      tagExcludeFileHandle.write(str(i)+"\n")
  tagExcludeFileHandle.close() #don't forget to close
  parsedBarcode=len(tagFractions.iloc[0])
  filteredBarcode=tcol
  print >>sys.stderr, "parsedBarcode=", parsedBarcode
  print >>sys.stderr, "filteredBarcode=", filteredBarcode
  print >>sys.stderr, "number of excluded >maxTagFraction tags:", parsedBarcode-filteredBarcode
  tagcol2tagidx = tagcol2tagidx[:tcol] #remove trailing -1 entries in tagcol2tagidx resulting from freq filter 
  #print "tagidx2tagcol",tagidx2tagcol
  frag2tagidx = list(fragdata.iloc[:,3]-1) #0-based fragment to tagidx lookup
  #print "frag2tagidx",frag2tagidx
  
  if debug:
    gtk1_time = time.time()
    print >>sys.stderr, "storing retro convertion from tcol to tag..."
    pandas.DataFrame(tagcol2tagidx, dtype=pandas.np.int64).to_csv(".".join([outputPrefix,"tagcol2tagidx","gz"])\
                                ,header=False, sep="\t", index=False, compression="gzip")
    gtk2_time = time.time()
    print >>sys.stderr, gtk2_time-gtk1_time, "s done"

  print >>sys.stderr, "creating matrix X..."
  #DONE: X(i,j) is overlaps of Region1(i) and Fragment with Tb(j)
  gtk1_time = time.time()

  X_i = list(overlap1.iloc[:,0]-1)              #grid: convert 1-based R index to to 0-based python index
  X_j = [tagidx2tagcol[frag2tagidx[i]] for i in overlap1.iloc[:,1]-1]
  X_data = pandas.DataFrame({ 'i' : X_i, 'j' : X_j }, dtype=pandas.np.int32) # region-i, tag-j
  X_data = X_data[ X_data['j']>0 ]
  
  if debug:
    print >>sys.stderr, "overlap1:", len(overlap1)
    print >>sys.stderr, "X_i:", len(X_i)
    print >>sys.stderr, "X_j:", len(X_j)
    print >>sys.stderr, "X_data (before skip cols):", len(X_data)
    print >>sys.stderr, "X_data (before group):", len(X_data)
    X_data.to_csv(".".join([outputPrefix,"X_df","gz"]), header=True, sep="\t", index=False, compression="gzip")

  X_data = X_data.groupby(['i','j']).size() #potentail problematic
  if debug: 
    print >>sys.stderr, "X_data (before reset index):", len(X_data)
  X_data = X_data.reset_index()
  if debug:
    print >>sys.stderr, "X_data (final):", len(X_data)
  X_data.columns = ['i', 'j', 'v']
  X_shape = (len(grid1), len(tagcol2tagidx))
  X = scipy.sparse.coo_matrix((X_data['v'], (X_data['i'], X_data['j'])), shape=X_shape)
  X = X.tocsr()
  print >>sys.stderr, "X.shape", X.shape
  print >>sys.stderr, "X.nnz", X.nnz
  if debug: 
    print >>sys.stderr, "saving matrix X..."
    X_data.to_csv(".".join([outputPrefix,"X_data","gz"]), sep="\t", index=False, compression="gzip")
    lib10x.save_sparse_csr(".".join([outputPrefix,"X"]),X)

  del X_i, X_j, X_data
  gc.collect()
  gtk2_time = time.time()
  print >>sys.stderr, gtk2_time-gtk1_time, "s done"
  
  print >>sys.stderr, "creating matrix Y..."
  #DONE: Y(j,k) is indicator if overlap of Region2(k) and Fragment with Tb(j)
  gtk1_time = time.time()
  if gridFile2 == gridFile1:
    Y_j = None
    Y_k = None
    Y_data = None
    Y = X.transpose()
    Y = Y.tocsr()
  else:
    Y_j = [tagidx2tagcol[frag2tagidx[i]] for i in overlap2.iloc[:,1]-1]
    Y_k = list(overlap2.iloc[:,0]-1)   #convert to 0-based python index
    Y_data = pandas.DataFrame({ 'j' : Y_j, 'k' : Y_k },dtype=pandas.np.int32) # region-i, tag-j
    Y_data = Y_data[ Y_data['j']>0 ]
    Y_data = Y_data.groupby(['j', 'k']).size()
    Y_data = Y_data.reset_index()
    Y_data.columns = ['j', 'k', 'v']
    Y_shape = (len(tagcol2tagidx), len(grid2))
    Y = scipy.sparse.coo_matrix((Y_data['v'], (Y_data['j'], Y_data['k'])), shape=Y_shape)
    Y = Y.tocsr()
  #print >>sys.stderr, len(Y_data)
  
  Y.data = Y.data>0
  print >>sys.stderr, "Y.shape", Y.shape
  print >>sys.stderr, "Y.nnz", Y.nnz
  if debug:
    print >>sys.stderr, "saving matrix Y..."
    if Y_data: Y_data.to_csv(".".join([outputPrefix,"Y_data","gz"]), sep="\t", index=False, compression="gzip")
    lib10x.save_sparse_csr(".".join([outputPrefix,"Y"]),Y)
    #scipy.io.mmwrite(".".join([outputPrefix,"Y"]), Y.toarray())
  del Y_j, Y_k, Y_data
  gc.collect()
  gtk2_time = time.time()
  print >>sys.stderr, gtk2_time-gtk1_time, "s done"
  
  #DONE create a sparse matrix M = X * t(Y)
  #NOTE: X and Y are big sparse matrices of 60K * 50K with around 30M entries each
  #      M is big dense matrix of 60K * 60K
  #      to save M, we probably can save X and Y and reconstruct M
  #      try save X and Y into .mtx first
  #      try save X and Y into .npy second 
  print >>sys.stderr, "controlling FDR..."
  #compoute tcut from FDR control: affordable computations N_c = M^2 * (1-ppois(u0,t))
  #cutFDR = 0
  #cutDUP = 0
  cutFDR = scipy.stats.poisson.ppf(1-expectFalseDiscover/float(X.shape[0]*Y.shape[1]),expectShareFrag)
  print >>sys.stderr, "cutFDR:", cutFDR
  cutDUP = scipy.stats.poisson.ppf(1-expectFalseDiscover/float(X.shape[0]*Y.shape[1]),expectJunctionFrag)
  print >>sys.stderr, "cutDUP:", cutDUP
  print >>sys.stderr, "multiplying matrix X and Y..."
  
  #DEBUG Block Dot Multiplication for Sparse Matrix
  #M1 = X.dot(Y)
  #print >>sys.stderr, "M1.shape", M1.shape
  #print >>sys.stderr, "M1.nnz:", M1.nnz
  #print >>sys.stderr, "M1.dtype:", M1.dtype
  #print >>sys.stderr, "M1.format:", M1.getformat()
  
  #M1.data[M1.data<=cutFDR]=0
  #M1.eliminate_zeros()
  #print >>sys.stderr, "M1.nnz after taking out <cutFDR:", M1.nnz
  
  sub_size_x = numpy.ceil(X.shape[0]/numSub)
  sub_size_y = numpy.ceil(Y.shape[1]/numSub)
  print >>sys.stderr, "sub_size_x", sub_size_x
  print >>sys.stderr, "sub_size_y", sub_size_y
  M = [ [ None for j in xrange(0,numSub) ] for i in xrange(0,numSub) ]
  #make sure the initialized lists are separate
  #print >>sys.stderr, "M=", M
  
  for i in xrange(0,numSub):
    for j in xrange(0,numSub):
      gtk1_time = time.time()
      print >>sys.stderr, "X-", i, ".Y-", j, ":",  
      M[i][j] = lib10x.sparse_sub_dot(X, Y, sub_size_x, i, sub_size_y, j)
      print >>sys.stderr, M[i][j].shape
      if debug: print >>sys.stderr, "nnzs:", M[i][j].nnz
      M[i][j].data[M[i][j].data<=cutFDR]=0
      M[i][j].eliminate_zeros()
      if debug: print >>sys.stderr, "nnz after taking out <cutFDR:", M[i][j].nnz
      M[i][j] = M[i][j].tocoo()
      gtk2_time = time.time()
      print >>sys.stderr, gtk2_time-gtk1_time, "s done"
  
  del X, Y, overlap1, overlap2, fragdata, frag2tagidx
  gc.collect()
  
  for i in xrange(0,numSub):
    M[i] = scipy.sparse.hstack(M[i])
  M = scipy.sparse.vstack(M)
  M = M.tocsr()

  #print >>sys.stderr, "M.format:", M.getformat()
  #print (M1!=M).nnz
  #print numpy.array_equal(M1.toarray(), M.toarray())
  #  #M = X.dot(Y)
  #  M_median_mad=lib10x.csr_median_mad(M)
  #  M_median=M_median_mad[0]
  #  M_mad=M_median_mad[1]
  #  #M_median = numpy.median(M.data.tolist()+[0]*(M.shape[0]*M.shape[1]))
  #  #M_mad = lib10x.mad(M.data.tolist()+[0]*(M.shape[0]*M.shape[1]))
  #  M_mean = M.mean()
  #  M_sd = lib10x.sparse_std(M)
  #  print >>sys.stderr, "M_mean:", M_mean
  #  print >>sys.stderr, "M_sd:", M_sd
  #  print >>sys.stderr, "M_median:", M_median
  #  print >>sys.stderr, "M_mad:", M_mad
  
  print >>sys.stderr, "M.shape:", M.shape
  print >>sys.stderr, "M.nnz:", M.nnz
  
  for row in skipIndex1:
    lib10x.csr_row_set_nz_to_val(M, row) 
  M.eliminate_zeros()
  print >>sys.stderr, "nnz after removing skip grid1:", M.nnz
  M = M.tocsc()
  for col in skipIndex2:
    lib10x.csc_col_set_nz_to_val(M, col) 
  M.eliminate_zeros()
  print >>sys.stderr, "nnz after removing skip grid2:", M.nnz
  
  if debug:
    print >>sys.stderr, "saving matrix M..."
    lib10x.save_sparse_csr(".".join([outputPrefix,"M"]),M)
  
  M = M.tocoo()
  print >>sys.stderr, "M is symmetric: ", lib10x.is_symmetric(M)
  M_data = pandas.DataFrame({ 'grid1' : M.row, 'grid2' : M.col, 'm12' : M.data }) #0-based grid indicies
  sel = M_data['grid1'] <= M_data['grid2']
  M_data = M_data.loc[sel,]
  M_data = M_data.sort_values(by=['m12','grid1','grid2'], ascending=[0,1,1]) #sorted
  M_data = M_data.reindex()
  
  #output possible duplications
  M_data_dup_idx = numpy.logical_and(numpy.abs(M_data['grid1']-M_data['grid2'])==0, M_data['m12']>cutDUP).as_matrix()
  #print M_data_dup_idx
  #print M_data_dup_idx.as_matrix()
  #print M_data_dup_idx.tolist()
  #print M_data_dup_idx.info()
  
  #M_data_dup_idx = M_data_dup_idx.reindex()
  #print M_data_dup_idx
  if debug: print >>sys.stderr, "maximum possible duplication grids: ", sum(M_data_dup_idx)
  c = 0
  for i in xrange(0,len(M_data_dup_idx)):
    if M_data_dup_idx[i]:
      c += 1
      if c > maxDuplication:
        M_data_dup_idx[i] = False
  if debug: print >>sys.stderr, "limiting duplication grids to: ", c
  M_data_dup = M_data[ M_data_dup_idx ]
  if debug: print >>sys.stderr, "writing results to %s ..." % ".".join([outputPrefix,"M_data","dup","gz"])
  M_data_dup.to_csv(".".join([outputPrefix,"dup","gz"]), sep="\t", index=False, compression="gzip")
  
  dfs = [grid1.iloc[M_data.iloc[M_data_dup_idx,0],:],grid2.iloc[M_data.iloc[M_data_dup_idx,1]],M_data.iloc[M_data_dup_idx,2].to_frame()]
  for d in dfs: d.reset_index(drop=True, inplace=True) #very important
  M_data_dup_bedpe=pandas.concat(dfs, names=bedpeLabels, axis=1)
  #for debug only
  #save = shelve.open("test.shelve")
  #save["grid1"] = grid1
  #save["grid2"] = grid2
  #save["M_data"] = M_data
  #save.close()
  #print type(grid1.iloc[M_data.iloc[M_data_dup_idx,0],:])
  #print type(grid2.iloc[M_data.iloc[M_data_dup_idx,0],:])
  #print type(M_data.iloc[M_data_dup_idx,2])
  #print type(M_data.iloc[M_data_dup_idx,2].to_frame())
  #M_data_dup_bedpe=pandas.concat([grid1.iloc[M_data.iloc[M_data_dup_idx,0],:],grid2.iloc[M_data.iloc[M_data_dup_idx,1]],M_data.iloc[M_data_dup_idx,2].to_frame()],names=bedpeLabels,axis=1) #use concat if and only if all are data frames
  #print(M_data.iloc[M_data_dup_idx,2])
  #M_data_dup_bedpe["m12"] = M_data.iloc[M_data_dup_idx,2]

  print >>sys.stderr, "writing results to %s ..." % ".".join([outputPrefix,"dup","bedpe"])
  M_data_dup_bedpe.to_csv(open(".".join([outputPrefix,"dup","bedpe"]),'w'), sep="\t", index=False, header=False)
  
  M_data_jct_idx = numpy.logical_and(M_data['grid2']-M_data['grid1']>3, M_data['m12']>cutFDR).as_matrix()
  #output possible junctions, only keep grid1<grid2 entries
  #filter out diagonal (+/-)3 elements
  #write M_data to csv
  #convert M_data to bed
  if debug: print >>sys.stderr, "maximum possible junction grids: ", sum(M_data_jct_idx)
  c = 0
  for i in xrange(0,len(M_data_jct_idx)):
    if M_data_jct_idx[i]:
      c += 1
      if c > maxJunction:
        M_data_jct_idx[i] = False
  if debug: print >>sys.stderr, "limiting junction grids to: ", c
  M_data_jct = M_data[ M_data_jct_idx ]
  if debug: print >>sys.stderr, "writing results to %s ..." % ".".join([outputPrefix,"jct","gz"])
  M_data_jct.to_csv(".".join([outputPrefix,"jct","gz"]), sep="\t", index=False, compression="gzip")
  
  dfs = [grid1.iloc[M_data.iloc[M_data_jct_idx,0],range(0,3)],grid2.iloc[M_data.iloc[M_data_jct_idx,1],range(0,3)],M_data.iloc[M_data_jct_idx,2].to_frame()]
  #print dfs
  for d in dfs: d.reset_index(drop=True, inplace=True) #very important
  M_data_dup_bedpe=pandas.concat(dfs, names=bedpeLabels, axis=1)
  print >>sys.stderr, "writing results to %s ..." % ".".join([outputPrefix,"jct","bedpe"])
  M_data_dup_bedpe.to_csv(open(".".join([outputPrefix,"jct","bedpe"]),'w'), sep="\t", index=False, header=False)
  
  cfgFile = "{}.{}".format(outputPrefix,"gridScan.cfg.json")
  scanConfig['G'] = genomeSize
  scanConfig['L'] = gridSize
  scanConfig['aNF'] = fragNumber           #molecule number regardless of barcode
  scanConfig['mNF'] = fragdataSize         #molecule number barcode > 2 molecules
  scanConfig['mLf'] = averageFragLength
  scanConfig['mNb'] = parsedBarcode
  scanConfig['fNb'] = filteredBarcode
  scanConfig['Tb'] = tagFractionSquareSum 
  scanConfig['Cem'] = fragNumber * averageFragLength / float(genomeSize)
  scanConfig['E0'] = parsedBarcode**2 * (gridSize/float(genomeSize))**2 * (tagFractionSquareSum)
  with open(cfgFile, 'w') as f: json.dump(scanConfig, f)
  
  if not debug:
    lib10x.rmTemporaryDir(tempdir)

#trouble shooting ValueError: negative dimension not allowed
#basically scipy version too low and int32 index overflow
#import lib10x
#X = lib10x.load_sparse_csr("/home/lixia/ws/ls411n/LS411N.wg10.bam.grid100K.test.X.npz")
#Y = lib10x.load_sparse_csr("/home/lixia/ws/ls411n/LS411N.wg10.bam.grid100K.test.Y.npz")

#DEBUG: 
#  test computation limitations
#  X = scipy.sparse.rand(60000, 60000, density=0.002, format='csr')
#  Y = X.transpose(copy=True)
#  M = X.dot(Y)   #memory bust

#inc_idx = [ x for x in xrange(0,len(a)) if a[x] not in exc ] #works
#print skipIndex1
#print len(M_data)
#print M_data['grid1'][0]
#print M_data['grid1'][len(M_data)-1]
#print [ x for x in xrange(0,len(M_data)) if M_data['grid1'].iloc[x] not in skipIndex1 ]
#quit()
#M_data = [ x for x in xrange(0,len(M_data)) if M_data['grid1'].iloc[x] not in skipIndex1 ] 
#M_data = [ x for x in xrange(0,len(M_data)) if M_data['grid2'].iloc[x] not in skipIndex2 ] 
#if True:
#
#  from numpy import mean, absolute
#  def mad(data, axis=None):
#    return numpy.median(absolute(data - numpy.median(data, axis)), axis)
#
#  #max_off = min(15, M.shape[0], M.shape[1])
#  max_off = 15
#
#  M_diagonal_median = [numpy.median( M_data[numpy.abs(M_data['grid1']-M_data['grid2']) == i]['m12'] ) for i in xrange(0,max_off)]
#  M_diagonal_mad = [mad( M_data[numpy.abs(M_data['grid1']-M_data['grid2']) == i]['m12'] ) for i in xrange(0,max_off)]
#
#  M_diagonal_mean = [numpy.mean( M_data[numpy.abs(M_data['grid1']-M_data['grid2']) == i]['m12'] ) for i in xrange(0,max_off)]
#  M_diagonal_std = [numpy.std( M_data[numpy.abs(M_data['grid1']-M_data['grid2']) == i]['m12'] ) for i in xrange(0,max_off)]
#
#  print >>sys.stderr, "M_diagonal_median:", M_diagonal_median
#  print >>sys.stderr, "M_diagonal_mad:", M_diagonal_mad
#  print >>sys.stderr, "M_diagonal_mean:", M_diagonal_mean
#  print >>sys.stderr, "M_diagonal_std:", M_diagonal_std
#if debug:
#
#  #max_off = min(15, M.shape[0], M.shape[1])
#  max_off = 15
#
#  #this is totally wrong. We cann't take median of M_data here
#  M_diagonal_median = [numpy.median( M_data[numpy.abs(M_data['grid1']-M_data['grid2']) == i]['m12'] ) for i in xrange(0,max_off)]
#  M_diagonal_mad = [mad( M_data[numpy.abs(M_data['grid1']-M_data['grid2']) == i]['m12'] ) for i in xrange(0,max_off)]
#
#  M_diagonal_mean = [numpy.mean( M_data[numpy.abs(M_data['grid1']-M_data['grid2']) == i]['m12'] ) for i in xrange(0,max_off)]
#  M_diagonal_std = [numpy.std( M_data[numpy.abs(M_data['grid1']-M_data['grid2']) == i]['m12'] ) for i in xrange(0,max_off)]
#
#  print >>sys.stderr, "M_diagonal_median:", M_diagonal_median
#  print >>sys.stderr, "M_diagonal_mad:", M_diagonal_mad
#  print >>sys.stderr, "M_diagonal_mean:", M_diagonal_mean
#  print >>sys.stderr, "M_diagonal_std:", M_diagonal_std
#
#  M_median = max(M_diagonal_median[3:])
#  M_mad = max(M_diagonal_mad[3:])
#  M_cutoff_norm_3mad = M_median + 3 * M_mad
#  M_cutoff_norm_FDR_empr = scipy.stats.norm.ppf(1-expectFalseDiscover/float(M.shape[0]**2),M_median,M_mad)
#  M_cutoff_pois_FDR_theo = scipy.stats.poisson.ppf(1-expectFalseDiscover/float(M.shape[0]**2),expectShareFrag)
#  M_cutoff_pois_FDR_empr = scipy.stats.poisson.ppf(1-expectFalseDiscover/float(M.shape[0]**2),M_median)
#
#  print >>sys.stderr, "M_cutoff_norm_3mad:", M_cutoff_norm_3mad
#  print >>sys.stderr, "M_cutoff_norm_FDR_empr:", M_cutoff_norm_FDR_empr
#  print >>sys.stderr, "M_cutoff_pois_FDR_theo:", M_cutoff_pois_FDR_theo
#  print >>sys.stderr, "M_cutoff_pois_FDR_empr:", M_cutoff_pois_FDR_empr
#
#  cutoff=max(M_cutoff_norm_3mad, M_cutoff_norm_FDR_empr, M_cutoff_pois_FDR_theo, M_cutoff_pois_FDR_empr)
#
#print >>sys.stderr, "cutoff:", cutoff
#print >>sys.stderr, "nnz after taking out <cutoff:", M.nnz
#M.data[M.data<=cutoff]=0
#M.eliminate_zeros()

#inc_idx = [ x for x in xrange(0,len(a)) if a[x] not in exc ] #works
#print skipIndex1
#print len(M_data)
#print M_data['grid1'][0]
#print M_data['grid1'][len(M_data)-1]
#print [ x for x in xrange(0,len(M_data)) if M_data['grid1'].iloc[x] not in skipIndex1 ]
#quit()
#M_data = [ x for x in xrange(0,len(M_data)) if M_data['grid1'].iloc[x] not in skipIndex1 ] 
#M_data = [ x for x in xrange(0,len(M_data)) if M_data['grid2'].iloc[x] not in skipIndex2 ] 

def run():
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("gridFile1",help="1st set of grids BED formatted ")
  parser.add_argument("fragFile",help="10x fragments BED formatted")
  parser.add_argument("hdf5File",help="10x statistics HDF5 formatted")
  parser.add_argument("outputPrefix",help="output prefix")
  parser.add_argument("-b","--gridFile2",default="same",help="2nd set of grids BED formatted (optional)")
  parser.add_argument("-s","--skipFile1",default="none",help="1st skipped grid regions BED formatted (optional)")
  parser.add_argument("-t","--skipFile2",default="same",help="2nd skipped grid regions BED formatted (optional)")
  parser.add_argument("-e","--expectFalseDiscover",type=int,default=10,help="expected false discovery allowed")
  parser.add_argument("-x","--maxDuplication",type=int,default=numpy.iinfo(numpy.int16).max,help="maximum expected number of >gridSize duplications")
  parser.add_argument("-y","--maxJunction",type=int,default=numpy.iinfo(numpy.int16).max,help="maximum expected number of >gridSize junctions")
  parser.add_argument("-p","--aneuploidyAllowance",type=float,default=4,help="adjust background rate by known aneuploidy")
  parser.add_argument("-n","--numSub",default=2,type=int,help="number of submatrix used for huge sparse matrices dot multiplication, memory=1/n^2")
  parser.add_argument("-m","--maxTagFraction",type=float,default=-1,help="max tag fraction to be considerrred")
  parser.add_argument("-g","--genomeSize",type=int,default=3200000000,help="bp size of the whole genome")
  parser.add_argument("-r","--gridSize",type=int,default=10000,help="bp size of the unit grid")
  parser.add_argument("-a","--debug",default=False,action="store_true",help="in debug mode")
  parser.add_argument('-d', '--tempdir', dest='tempdir',\
    default=os.path.join(os.environ['HOME'],'zoomx_'+lib10x.randascii(4)),\
    help="root dir for keeping temporary files, default (last 4 digits random): %(default)s")
  try:
    args = parser.parse_args()
  except SystemExit:
    args = None
  main(args)
  print >>sys.stderr, "GRIDSCANDONE; elapsed time", time.time()-start_time, "seconds"

if __name__ == '__main__':
  run()
