#!/usr/bin/env python
#this script uses two spaces as indent
import os, sys, time, gc, argparse, re, tempfile, string
import pysam
import numpy, pandas, json, random, shelve, sets
import scipy.sparse
import scipy.sparse.csgraph

try:
  import lib10x
  print >>sys.stderr, "loaded downloaded zoomx"
except ImportError:
  import zoomx.lib10x
  print >>sys.stderr, "load installed zoomx"
start_time=time.time()

def main(args):

  """
    DESIGN:
      Call:
        (chr1, st1, ed1, ps1, chr2, st2, ed2, ps2, CallType, CallID, SupportIDs)
      Support:
        (chr1, st1, ed1, dir1, chr2, st2, ed2, dir2, CallID, DataString)
      Data:                               DataString:
        (PairType, PairName, HapBlock) =>   PairType;PairName;HapBlock
      PairType:
        0 = Normal, 1 = Insert, 2 = Split, ...
      
    INPUT:
      bam
      bam.cov
      bam.M_data
      call.bedpe
      impp.bedpe
      disc.bedpe
      splt.bedpe
        
    OUTPUT:
      call.bedpe (merged, with type, id, genotype)
      call.pdf (mol-plots and heat-plots) 
      
  """
  
  if args:
    outPrefix = args.outPrefix
    bamFile = args.bam
    callFile = args.call
    cfgFile = args.cfg
    gapFile = args.gapFile
    frBedpeFile = args.frBedpeFile
    rfBedpeFile = args.rfBedpeFile
    ffBedpeFile = args.ffBedpeFile
    rrBedpeFile = args.rrBedpeFile 
    moleBedFile = args.moleBedFile 
    debug = args.debug
    mergeOnly = args.mergeOnly
    tempdir = args.tempdir
    numprocs = args.numprocs
    altExpect = args.altExpect
    forcedExpect = args.forcedExpect
    maxNum = args.maxNum
  else:
    bamFile = "LS411N_CR1.wg12.bam"
    if not os.path.isfile(bamFile): quit()
    print >>sys.stderr, "fall back to run program tests..."
    outPrefix = "none"
    callFile = "LS411N_CR1.wg12.bam.jct.bedpe"
    cfgFile = "LS411N_CR1.wg12.bam.gridScan.cfg.json"
    gapFile = "none"
    frBedpeFile = "LS411N_CR1.wg12.bam.FR.bedpe"
    rfBedpeFile = "LS411N_CR1.wg12.bam.RF.bedpe"
    ffBedpeFile = "LS411N_CR1.wg12.bam.FF.bedpe"
    rrBedpeFile = "LS411N_CR1.wg12.bam.RR.bedpe"
    moleBedFile = "LS411N_CR1.wg12.bam.multiFragBar.bed.a0"
    tempdir = 'zoomx_'+lib10x.randascii(4)
    numprocs = 2
    debug = False
    mergeOnly = False
    altExpect = 25
    forcedExpect = 10
    maxNum = 10
  
  lib10x.setTemporaryDir(tempdir)
  with open(cfgFile, 'r') as f:
    scanConfig = json.load(f)
  minMolCnt = numpy.ceil((1-altExpect/100.)*scanConfig['E0'] + altExpect/100.*scanConfig['Cem'])
  if forcedExpect != 0:
    minMolCnt = max(minMolCnt, forcedExpect)  
  sys.stderr.write("minMolCnt="+str(minMolCnt)+"\n")
  
  if debug:
    shelve_file = os.path.join(tempfile.tempdir,'shelve')
    sys.stderr.write("shelve_file='"+str(shelve_file)+"'\n")
    save = shelve.open(shelve_file) #shelve file 
  
  if outPrefix == "none":
    outPrefix = callFile
  if frBedpeFile == "none":
    frBedpeFile = ".".join([bamFile,'FR','bedpe'])
  if rfBedpeFile == "none":
    rfBedpeFile = ".".join([bamFile,'RF','bedpe'])
  if ffBedpeFile == "none":
    ffBedpeFile = ".".join([bamFile,'FF','bedpe'])
  if rrBedpeFile == "none":
    rrBedpeFile = ".".join([bamFile,'RR','bedpe'])
  if moleBedFile == "none":
    moleBedFile = ".".join([bamFile,'multiFragBar','bed','a0'])
  
  sys.stderr.write("refineBreak.py -m {} -f {} -e {} -n {} {} {} {}\n".format(maxNum,forcedExpect,altExpect,numprocs,bamFile,callFile,cfgFile)) 
  sys.stderr.write("args="+str(args)+"\n")
  if debug:
    sys.stderr.write("outPrefix="+str(outPrefix)+"\n")
    sys.stderr.write("bamFile="+str(bamFile)+"\n")
    sys.stderr.write("callFile="+str(callFile)+"\n")
    sys.stderr.write("frBedpeFile="+str(frBedpeFile)+"\n")
    sys.stderr.write("rfBedpeFile="+str(rfBedpeFile)+"\n")
    sys.stderr.write("ffBedpeFile="+str(ffBedpeFile)+"\n")
    sys.stderr.write("rrBedpeFile="+str(rrBedpeFile)+"\n")
    sys.stderr.write("moleBedFile="+str(moleBedFile)+"\n")
    sys.stderr.write("debug="+str(debug)+"\n")
    sys.stderr.write("mergeOnly="+str(mergeOnly)+"\n")
    sys.stderr.write("tempdir="+str(tempfile.tempdir)+"\n")
  
  callLabels=["chr1","start1","end1","chr2","start2","end2","barcnt"]
  callTypes=["string",numpy.int64,numpy.int64,"string",numpy.int64,numpy.int64,numpy.int64]
  mergeLabels=["chr1","start1","end1","chr2","start2","end2","name","barcnt"]
  mergeTypes=["string",numpy.int64,numpy.int64,"string",numpy.int64,numpy.int64,"string",numpy.int64]
  rdpLabels=["chr1","start1","end1","chr2","start2","end2","hp"]
  rdpTypes=["string",numpy.int64,numpy.int64,"string",numpy.int64,numpy.int64,numpy.int64]
  bed1Labels=["chr1","start1","end1","barcnt","index"]
  bed2Labels=["chr2","start2","end2","barcnt","index"]
  clust1Labels=["chr1","start1","end1","barcnt","index","clust"]
  clust2Labels=["chr2","start2","end2","barcnt","index","clust"]
  bedLabels=["chr","start","end","barcnt"]
  
  #rawCalls=pandas.read_table(callFile, names=callLabels)
  rawCalls=pandas.read_table(callFile,names=callLabels,dtype=dict(zip(callLabels,callTypes))) #only read cols 1-7
  rawCalls=rawCalls[rawCalls['barcnt']>minMolCnt] #keep only min barcnt and above
  if (len(rawCalls) > maxNum): rawCalls=rawCalls[0:maxNum]
  sys.stderr.write("maxNum={}\n".format(maxNum)) #reduce rawCalls to maxNum allowed
  sys.stderr.write("rawCalls={}\n".format(len(rawCalls)))

  tmpBedpeFile=None
  if (len(rawCalls) > 0):
    #save tmpBedpefile for pairtopair overlapping
    tmpRawFile=lib10x.getTemporaryFile('bedpe_')
    if debug: sys.stderr.write("tmpRawFile='"+str(tmpRawFile)+"'\n")
    rawCalls.to_csv(tmpRawFile,header=False, sep="\t", index=False)

    if (gapFile != "none"): #this removed additional rawCalls
      sys.stderr.write("filtering gap calls...\n")
      #do gapFile-based filtering of calls
      #assuming gapFile has already been padded
      tmpBedpeFile=lib10x.getTemporaryFile('bedpe_')
      cutString=",".join([str(x) for x in xrange(1,rawCalls.shape[1]+1)]) #include the last index column
      gapfilter_cmd = "bedtools pairtobed -type neither -a %s -b %s | cut -f%s >%s" % (tmpRawFile, gapFile, cutString, tmpBedpeFile)
      if lib10x.run_cmd(gapfilter_cmd) == False:
        sys.stderr.write("Error in running gapfilter overlapping..."); quit()
    else:
      tmpBedpeFile=tmpRawFile

    rawCalls=pandas.read_table(tmpBedpeFile, names=callLabels,dtype=dict(zip(callLabels,callTypes))) #reload after filter
    rawCalls['index']=xrange(0,rawCalls.shape[0]) #add index to bedpe file, index is fixed here
    rawCalls.to_csv(tmpBedpeFile,header=False, sep="\t", index=False) #no more add/remove entries of rawCalls from here
    if debug: sys.stderr.write("tmpBedpeFile='"+str(tmpBedpeFile)+"'\n")

  mergeBedpe=pandas.DataFrame(dict(zip(mergeLabels,[ [] for i in range(len(mergeLabels)) ])), columns=mergeLabels)
  if (len(rawCalls) > 0):
    nRawCalls = len(rawCalls)
    #pairtopair overlapping
    sys.stderr.write("merging raw calls...\n")
    sys.stderr.write("nRawCalls={}\n".format(nRawCalls))
    tmpOverlapFile=lib10x.getTemporaryFile('bedpe_')
    pair2pair_cmd="bedtools pairtopair -type both -is -slop 5 -a %s -b %s | cut -f8,16 >%s" % (tmpBedpeFile, tmpBedpeFile, tmpOverlapFile)
    if debug: sys.stderr.write("pair2pair_cmd="+str(pair2pair_cmd)+"\n")
    if debug: sys.stderr.write("tmpOverlapFile='"+str(tmpOverlapFile)+"'\n")
    if lib10x.run_cmd(pair2pair_cmd) == False:
      sys.stderr.write("Error in running pair2pair overlapping..."); quit()

    end_time=time.time()
    sys.stderr.write("overlapping done:{}s\n".format(end_time-start_time))

    sys.stderr.write("loading overlapping graph...")
    overlapCalls=pandas.read_table(tmpOverlapFile, names=["bedpe1","bedpe2"], dtype={'bedpe1':numpy.int64,'bedpe2':numpy.int64})
    #load pairtopair overlapping into a sparse matrix
    if debug:
      save['overlapCalls']=overlapCalls
    ones=[1]*overlapCalls.shape[0]
    sys.stderr.write("analyzing connecting component...")
    overlapGraph=scipy.sparse.coo_matrix((ones,(overlapCalls.iloc[:,0],overlapCalls.iloc[:,1])), shape=(nRawCalls,nRawCalls))
    if debug:
      save['overlapGraph']=overlapGraph
    #find connected components; only keep correctly ordered components
    numComponents, components = scipy.sparse.csgraph.connected_components(overlapGraph)
    #find connected entries of rawCalls
    #c -> (c/23, c%23)
    #reduce connected components to a new bedpe format

    mergeChr1=numpy.array([None] * numComponents, dtype="object")
    mergeStart1=numpy.array([None] * numComponents, dtype="object")
    mergeEnd1=numpy.array([None] * numComponents, dtype="object")
    mergeChr2=numpy.array([None] * numComponents, dtype="object")
    mergeStart2=numpy.array([None] * numComponents, dtype="object")
    mergeEnd2=numpy.array([None] * numComponents, dtype="object")
    mergeMolCnt=numpy.array([None] * numComponents, dtype="object")
    for ci in range(numComponents): # you might get flipped 
      ri=numpy.where(components == ci)
      mergeChr1[ci]=list(numpy.unique(rawCalls.iloc[ri]['chr1']))[0]
      mergeStart1[ci]=min(rawCalls.iloc[ri]['start1'])
      mergeEnd1[ci]=max(rawCalls.iloc[ri]['end1'])
      mergeChr2[ci]=list(numpy.unique(rawCalls.iloc[ri]['chr2']))[0]
      mergeStart2[ci]=min(rawCalls.iloc[ri]['start2'])
      mergeEnd2[ci]=max(rawCalls.iloc[ri]['end2'])
      mergeMolCnt[ci]=max(rawCalls.iloc[ri]['barcnt']) #

    mergeName=["call_"+str(x) for x in xrange(1,len(mergeChr1)+1)]
    mergeBedpe=pandas.DataFrame(dict(zip(mergeLabels,[mergeChr1,mergeStart1,mergeEnd1,mergeChr2,mergeStart2,mergeEnd2,mergeName,mergeMolCnt])),\
      columns=mergeLabels) # natually sorted, good!
    mergeBedpe=mergeBedpe.sort_values(by=['barcnt'],ascending=False)

  mergeBedpeFile=".".join([outPrefix,"merged"])
  sys.stderr.write("output merged calls\n")
  mergeBedpe.to_csv(mergeBedpeFile,header=True, sep="\t", index=False)
  sys.stderr.write("mergedBedpeFile='"+str(mergeBedpeFile)+"'\n")
  end_time=time.time()
  sys.stderr.write("output merged done:{}s\n".format(end_time-start_time))
  
  if mergeOnly:   #quit now if only merging events, mergeBedpe has 8 cols
    quit()
  
  sys.stderr.write("attaching short read support\n")
  sys.stderr.write("merging FR rdp\n")
  frSupport = lib10x.ParallelIdentifySupportRDP(mergeBedpe, frBedpeFile, rdpLabels, rdpTypes, numprocs, debug)
  end_time=time.time()
  sys.stderr.write("merging FR done:{}s\n".format(end_time-start_time))

  sys.stderr.write("merging RF rdp\n")
  rfSupport = lib10x.ParallelIdentifySupportRDP(mergeBedpe, rfBedpeFile, rdpLabels, rdpTypes, numprocs, debug)
  end_time=time.time()
  sys.stderr.write("merging RF done:{}s\n".format(end_time-start_time))
  
  sys.stderr.write("merging FF rdp\n")
  ffSupport = lib10x.ParallelIdentifySupportRDP(mergeBedpe, ffBedpeFile, rdpLabels, rdpTypes, numprocs, debug)
  end_time=time.time()
  sys.stderr.write("merging FF done:{}s\n".format(end_time-start_time))
  
  sys.stderr.write("merging RR rdp\n")
  rrSupport = lib10x.ParallelIdentifySupportRDP(mergeBedpe, rrBedpeFile, rdpLabels, rdpTypes, numprocs, debug)
  end_time=time.time()
  sys.stderr.write("merging RR rdp:{}s\n".format(end_time-start_time))

  sys.stderr.write("attaching unique molecule support\n")
  sys.stderr.write("merging non-self-overlapping molecules\n")
  miSupport = lib10x.ParallelIdentifyNonSelfOverlapMole(mergeBedpe, moleBedFile, callLabels, callTypes, numprocs, debug)
  end_time=time.time()
  sys.stderr.write("merging done:{}s\n".format(end_time-start_time))
  
  mergeBedpe['fr']=frSupport
  mergeBedpe['rf']=rfSupport
  mergeBedpe['ff']=ffSupport
  mergeBedpe['rr']=rrSupport
  mergeBedpe['mi']=miSupport
  
  if debug:
    print(frSupport)
    print(rfSupport)
    print(ffSupport)
    print(rrSupport)
    print(miSupport)
    print(mergeBedpe)
  
  sys.stderr.write("output refined calls\n")
  refineBedpeFile=".".join([outPrefix,"refined"])
  #refineBedpe=mergeBedpe.sort_values(by=['barcnt'],ascending=False) #this is already sorted
  mergeBedpe.to_csv(refineBedpeFile,header=True, sep="\t", index=False)
  sys.stderr.write("refineBedpeFile='"+str(refineBedpeFile)+"'\n")

  if not debug:
    lib10x.rmTemporaryDir(tempdir)

#below obsolete

# find clusters and grids supporting them
#rawSort1=pybedtools.BedTool().from_dataframe(rawCalls.loc[:,bed1Labels]).sort().to_dataframe()
#rawSort2=pybedtools.BedTool().from_dataframe(rawCalls.loc[:,bed2Labels]).sort().to_dataframe()
#rawSort1.columns=bed1Labels
#rawSort2.columns=bed2Labels
#rawClust1=pybedtools.BedTool().from_dataframe(rawSort1).cluster().to_dataframe()
#rawClust2=pybedtools.BedTool().from_dataframe(rawSort2).cluster().to_dataframe()
#rawClust1.columns=clust1Labels
#rawClust2.columns=clust2Labels
#sortClust1=rawClust1.sort_values(by=['index'])
#sortClust2=rawClust2.sort_values(by=['index'])
#Clust1=list(set(sortClust1['clust'])) #left clusters
#Clust2=list(set(sortClust2['clust'])) #right clusters
#assert Clust1[len(Clust1)-1]==len(Clust1)
#assert Clust2[len(Clust2)-1]==len(Clust2)
#Clust=numpy.zeros(shape=(len(Clust1),len(Clust2)))
#ClustToRaw1=dict()
#ClustToRaw2=dict()
#for i in xrange(0, sortClust1.shape[0]):
#  assert sortClust1['index'].iloc[i] == i
#  k=sortClust1['clust'].iloc[i]-1
#  ClustToRaw1[k]=ClustToRaw1.get(k,[]) + [i]
#for i in xrange(0, sortClust2.shape[0]):
#  assert sortClust2['index'].iloc[i] == i
#  k=sortClust2['clust'].iloc[i]-1
#  ClustToRaw2[k]=ClustToRaw2.get(k,[]) + [i]
#NumClust1=Clust.shape[0]
#NumClust2=Clust.shape[1]
#NumClust=NumClust1+NumClust2
#Diag1=scipy.sparse.eye(NumClust1,dtype=numpy.int64,format='dok')
#Diag2=scipy.sparse.eye(NumClust2,dtype=numpy.int64,format='dok')
#ClustDOK=scipy.sparse.dok_matrix(Clust>=1,dtype=numpy.int64)
#GraphCOO=scipy.sparse.hstack([scipy.sparse.vstack([Diag1,ClustDOK]), scipy.sparse.vstack([ClustDOK.transpose(),Diag2])]) #this action converts DOK to COO
#NumGraphComponent,GraphComponent=scipy.sparse.csgraph.connected_components(GraphCOO)
## output clusters in bedpe format
## each component is an event
#EventComponent1=dict()
#EventComponent2=dict()
#for i in xrange(0,len(GraphComponent)): #Component of Graph to 
#  if i < NumClust1:
#    c = GraphComponent[i]
#    EventComponent1[c] = EventComponent1.get(c,set()).union(set(ClustToRaw1[c]))
#  else:
#    c = GraphComponent[i-NumClust1]
#    EventComponent2[c] = EventComponent2.get(c,set()).union(set(ClustToRaw2[c]))

# DFS component searching

#save['Clust']=Clust
#save['ClustToRaw1']=ClustToRaw1
#save['ClustToRaw2']=ClustToRaw2
#save['rawCalls']=rawCalls
#save['GraphCOO']=GraphCOO
#save['EventComponent1']=EventComponent1
#save['EventComponent2']=EventComponent2
#
#save.close()
#quit()
#
#rawChange1=numpy.where(rawClust1.values[range(1,rawClust1.shape[0]),4] - rawClust1.values[range(0,rawClust1.shape[0]-1),4])[0]
#rawChange2=numpy.where(rawClust2.values[range(1,rawClust2.shape[0]),4] - rawClust2.values[range(0,rawClust2.shape[0]-1),4])[0]
#RIGHT? real change point are points that in both rawChange1 and rawChange2
#rawChange=sorted(list(sets.Set(rawChange1.tolist()).intersection(sets.Set(rawChange2.tolist())))) # has to be sorted
#rawChangeStart=[0] + rawChange
#rawChangeEnd=[ x - 1 for x in rawChange ] + [rawClust1.shape[0]-1]
##conbine rawSort1, rawSort2 to mergeSort
#mergeChr1=rawSort1.iloc[rawChangeStart,0].values
#mergeStart1=rawSort1.iloc[rawChangeStart,1].values
#mergeEnd1=rawSort1.iloc[rawChangeEnd,2].values
#mergeChr2=rawSort2.iloc[rawChangeStart,0].values
#mergeStart2=rawSort2.iloc[rawChangeStart,1].values
#mergeEnd2=rawSort2.iloc[rawChangeEnd,2].values
#mergeMol=[ rawSort2.iloc[rawChangeStart[x]:(rawChangeEnd[x]+1),3] for x in range(0,len(rawChangeStart)) ]
#mergeMolCnt=[ numpy.mean(x.values) for x in mergeMol ]
#mergeSort=pandas.DataFrame(dict(zip(callLabels,[mergeChr1,mergeStart1,mergeEnd1,mergeChr2,mergeStart2,mergeEnd2,mergeMolCnt])),\
#  columns=callLabels)
#
#save['mergeSort']=mergeSort

#Now to add in more supporting read pairs to mergeSort

#Now output final result

def run():
  parser = argparse.ArgumentParser()
  parser.add_argument("bam", help="sorted 10x bamfile, in Longranger BAM format")
  parser.add_argument("call", help="SV call bedpefile, in bedpe format") 
  parser.add_argument("cfg",help="parameter config file, in JSON format")
  parser.add_argument("-g","--gapFile", default="none",help="bedfile of boundaries of known gaps >5K + slop(if any)")
  parser.add_argument("-m","--maxNum", default=numpy.iinfo(numpy.int64).max,type=int,
    help="max number of events to process (reduce this if time constrained)")
  parser.add_argument("-e","--altExpect", default=25,type=int, help="minimum mixture fraction for alternative, specify a percentage between 10-100, default 25")
  parser.add_argument("-f","--forcedExpect",default=10,type=int, help="forced minimum molecule support for SV events: default %(default)s")
  parser.add_argument("-i","--frBedpeFile", default="none",help="impp FR read file in bedpe, default=bam.FR.bedpe")
  parser.add_argument("-j","--ffBedpeFile", default="none",help="disc FF read file in bedpe, default=bam.FF.bedpe")
  parser.add_argument("-k","--rrBedpeFile", default="none",help="splt RR read file in bedpe, default=bam.RR.bedpe")
  parser.add_argument("-l","--rfBedpeFile", default="none",help="impp RF read file in bedpe, default=bam.RF.bedpe")
  parser.add_argument("-x","--moleBedFile", default="none",help="linked-read single molecule file, default=bam.multiFragBar.bed.a0")
  parser.add_argument("-o","--outPrefix",help="outputPrefix",default="none")
  parser.add_argument("-a","--debug",default=False,action="store_true",help="in debug mode")
  parser.add_argument("-b","--mergeOnly",default=False,action="store_true",help="in mergeOnly mode")
  parser.add_argument("-n","--numprocs",default=1,type=int,help="number of parallele process")
  parser.add_argument("-d","--tempdir", dest="tempdir",\
    default=os.path.join(os.environ['HOME'],'zoomx_'+lib10x.randascii(4)),\
    help="root dir for keeping temporary files, default (last 4 digits random): default: %(default)s")
  try:
    args = parser.parse_args()
  except SystemExit:
    args = None
  main(args)
  print >>sys.stderr, "REFINEBREAKDONE; elapsed time", time.time()-start_time, "seconds"

if __name__ == '__main__':
  run()
