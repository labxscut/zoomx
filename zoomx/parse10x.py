#!/usr/bin/env python
#using two spaces as indent

import time
import os
import pysam
import collections
import numpy
import bisect
import json
import pandas
import re
import sys
import argparse
import copy
import gc
from itertools import chain
#properly import zoomx modules
try:
  import lib10x
except ImportError:
  import zoomx.lib10x
start_time=time.time()

def main(args):

  """
		Parse10x Module of ZoomX for parsing 10X Genomics Bam
    
    INPUT:
      Bam:
        $pf.bam

    VARIABLE:
      Fragment:
        (gip, barcode, [fr1_gip, fr2_gip, ...]) 
      barfragdata:
        barfragdata[alpha]=dict(b:[fb1,fb2,...])
      fragdata:
        fragdata[alpha]=[f1,f2,...]
      multiFragBar:  
        multiFragBar i.e. barcodes having multiple molecules
        multiFragBar[alpha][i]=i-th barcode having len(barfragdata[barcode])>1 i.e. i-th barcode have more than one molecule
      fragCountMultiFragBar:
        fragCountMultiFragBar i.e. molecule counts for barcodes having multiple molecules
        fragCountMultiFragBar[alpha][i]=number of molecules having barcode multiFragBar[i]
      S(s):
        fragBarFragRdpCount[alpha][s]=number of molecules having a size (number of co-molecule rdps) s
      D(d):
        multiRdpFragRdpDistance[alpha][d]=number of d-spaced consecutive reads on a molecule
      T(t):
        fragFracMultiFragBar[alpha][i]=fraction of fragCountMultiFragBar[i] in total fragCountMultiFragBar
    
    OUTPUT:
      $pf.h5.aK:
        containing S(s), D(d) and T(t)
      $pf.[fragBar|multiFragBar]Bedfile.bed.aK:
        containing all molecules | mb-moleculess
        frag_chr  frag_begin  frag_end  frag_bar  frag_rdp_count  bar_frag_count
      
  """
  
  if args:
    outPrefix = args.outPrefix
    longRangerVersion = args.longRangerVersion 
    bamfile = args.bam
    covFile = args.covFile
    maxRead = args.maxRead
    setAlpha = [int(x) for x in args.setAlpha.split(',')]
    debug = args.debug
  else:
    bamfile="LS411N_CR1.wg12.bam"
    if not os.path.isfile(bamfile): quit() # quit if not in test directory
    print >>sys.stderr, "fall back to run program tests..."
    longRangerVersion = 1.2
    outPrefix = bamfile
    covFile="none"
    setAlpha=[0]
    maxRead=numpy.inf
    debug = False
    #maxRead=200000

  sys.stderr.write("parse10x.py -b {} {}\n".format(maxRead,bamfile))
  sys.stderr.write("args="+str(args)+"\n")
  #print >>sys.stderr, "maxRead=", maxRead
  if outPrefix == "none":
    outPrefix = bamfile
  
  imBedpeFile = open(outPrefix+".IM.bedpe",'w') 
  ffBedpeFile = open(outPrefix+".FF.bedpe",'w') 
  rrBedpeFile = open(outPrefix+".RR.bedpe",'w') 
  
  samfile = pysam.Samfile(bamfile, "rb")
  
  #1. read in the head file and construct the g2cp and cp2g maps
  refdata=samfile.header['SQ']
  # cp2g: g = A[chr]+pos; cp2g('chr',pos,A)
  # g2cp: chr = max_c : A[c]<g; pos=g-A[chr]
  chrname=map(lambda t : t['SN'], refdata)
  chrlen=map(lambda t : t['LN'], refdata)
  acclen=numpy.cumsum(chrlen)                                  
  A=collections.OrderedDict(zip(chrname,list(acclen)))                  #A: accumulated length dict (including the current chr)
  acclen_A=list(acclen)                                                 #acclen_A: accumulated length vec working with A
  C=collections.OrderedDict(zip(chrname,[1]+list(acclen[:-1]+1)))       #C: accumulated length dict (excluding the current chr)
  acclen_C=[1]+list(acclen[:-1]+1)                                      #acclen_C: accumulated length vec working with C
  refdata={'chrname':chrname,'chrlen':chrlen,'acclen_A':acclen_A,'acclen_C':acclen_C,'A':A,'C':C}
  if debug:
    print >>sys.stderr, "A=", A
    print >>sys.stderr, "acclen_A=", acclen_A
    print >>sys.stderr, "C=", C
    print >>sys.stderr, "acclen_C=", acclen_C
  
  #1. calculate the coverage 
  covVec=numpy.array([0]*int(acclen_A[-1]/100),dtype='float') 
  print >>sys.stderr, "allocated covVec:", len(covVec)
  cov_cut=numpy.inf
  cov_cut2=numpy.inf
  if covFile != "none":
    covLabels=['chr','start','end','counts']
    covTypes=dict(zip(covLabels,[str,numpy.int32,numpy.int32,numpy.int32]))
    covTable=pandas.read_csv(covFile,header=None,names=covLabels,dtype=covTypes,delimiter='\t')
    covChrChange=numpy.where(~covTable['chr'].duplicated())[0].tolist()
    print >>sys.stderr, "covChrChange:", covChrChange
    for i in xrange(0,len(covChrChange)):
      idxEnd = covChrChange[i+1] if i < len(covChrChange)-1 else len(covTable)
      idxStart = covChrChange[i]
      covStart=int(lib10x.cp2g(str(covTable['chr'][idxStart]),covTable['start'][idxStart],C)/100)
      covEnd=int(lib10x.cp2g(str(covTable['chr'][idxEnd-1]),covTable['start'][idxEnd-1],C)/100)
      #NOTE: coverage is read from covTable (.cov100) file and convert to global index position (gip/100) table
      #      so the coverage should't be wrong because of reference difference
      covVec[covStart:(covEnd+1)] = covTable['counts'][idxStart:idxEnd]
      #NOTE: covVec is an gip/100 indexed table
    
    cov_mean=numpy.mean(covVec[covVec!=0])
    cov_sd=numpy.std(covVec[covVec!=0])
    cov_cut2=cov_mean + 2*cov_sd                      #5 sd from mean, potentially type-1 error 2.3% 
    cov_median=numpy.median(covVec[covVec!=0])        
    cov_mad=lib10x.mad(covVec[covVec!=0])
    cov_cut=cov_median + 2*1.48*cov_mad
    print >>sys.stderr, "cov_median=", cov_median
    print >>sys.stderr, "cov_mad=", cov_mad
    print >>sys.stderr, "cov_mean=", cov_mean
    print >>sys.stderr, "cov_sd=", cov_sd
  
  print >>sys.stderr, "cov_cut=", cov_cut
  print >>sys.stderr, "cov_cut2=", cov_cut2
  
  #
  
  #2. define a python dict structure to cluster and hold the { b: (f11_idx, f12_idx, ...) } structure
  #3. define a python vector structure to hold [ (g(s11),b1,g(s11),g(s12),g(s1M)) ]
  #  3.1 molecule data [frag_start, barcode, (frag_rd1, ..., frag_rdM)]
  #  3.2 barcode data { bar : (frag_1, ..., frag_X) }
  barrpdata = dict()
  statReadMax = 1000000 #max number of read to load for insert and read len statistics
  rpinsert = []
  readlen = []
  N_all = 0
  N_qc = 0
  N_jc = 0
  N_b = 0
  N_goodLeftRead = 0
  
  if longRangerVersion < 1.2: #use nearest neighbor clustering
  
    print >>sys.stderr, "using manual nearest neighbour clustering (distance<alpha) to infer molecule..."
  
    numAlpha = len(setAlpha)
    barfragdata = dict(zip(setAlpha,[ dict() for x in xrange(numAlpha)])) #NOTE: * doesn't create independent dictionaries 
    fragdata = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    N_f =  dict(zip(setAlpha,[0]*numAlpha)) #alpha is alpha, here * does create independent integers
  
    for rd in samfile.fetch():
    
      if N_all % 100000 == 0 and N_all != 0: print >>sys.stderr, N_all, "reads processed"
    
      if N_all >= maxRead: break
    
      N_all += 1
    
      #filter based on if bearing barcode, throw Barcode unresolved read pairs
      if not dict(rd.tags).get('BX'):
        continue
  
      #filter based on sam flags
        # i.e. -f 0x001 -F 0x004 -F 0x008 -F 0x100 -F 0x800 -F 0x400 -F 0x200
        #throw unmapped and low quality read and qc failed read read pairs
      if (not rd.is_paired) or rd.is_unmapped or rd.mate_is_unmapped or rd.is_qcfail or rd.is_secondary or rd.is_supplementary or rd.is_duplicate:
        continue
  
      try:
         hp = int(dict(rd.tags)['HP'])
         ps = int(dict(rd.tags)['PS'])
      except:
         hp = -1
         ps = -1
    
      #output imBedpe/ffBedpe/rrBedpe
      #TODO: output split into ff, rr, im(1st=+,2nd=-) and mi(
      if ( not rd.is_reverse ) and rd.mate_is_reverse and ( not rd.is_proper_pair ) : #IM
        imBedpeFile.write('\t'.join([rd.reference_name, str(rd.reference_start), str(rd.reference_start+100), \
                         rd.next_reference_name, str(rd.next_reference_start), str(rd.next_reference_start+100), str(hp)])+'\n' )
  
      if ( not rd.is_reverse ) and ( not rd.mate_is_reverse ) : #FF
        ffBedpeFile.write('\t'.join([rd.reference_name, str(rd.reference_start), str(rd.reference_start+100), \
                         rd.next_reference_name, str(rd.next_reference_start), str(rd.next_reference_start+100), str(hp)])+'\n' )
  
      if rd.is_reverse and rd.mate_is_reverse: #RR
        rrBedpeFile.write('\t'.join([rd.reference_name, str(rd.reference_start), str(rd.reference_start+100), \
                         rd.next_reference_name, str(rd.next_reference_start), str(rd.next_reference_start+100), str(hp)])+'\n' )
  
      #filter based on read quality
      if rd.mapq<30:
        continue
    
      N_qc += 1
    
      if not rd.reference_id==rd.next_reference_id: #we only consider read pairs mapped to the same chromosome
        continue
    
      if rd.pos>rd.mpos: #we only consider the leftmost read of any read pair
        continue
    
      if rd.pos==rd.mpos and rd.is_read2: #if it cann't be distinguished by pos we distinguish by first/second read
        continue #NOTE: this make sure one read pair is only counted once and reduce excessive zeros in spacing
      #NOTE: pysam pos is 0-based !!!!
    
      if rd.pos<0 or rd.pos>=A[samfile.get_reference_name(rd.reference_id)]:
        continue #NOTE: if exceeding the bounds of reference length ignore
    
      bar = dict(rd.tags)['BX'] #BX is the barcode field
      gip_GoodReadLeft = lib10x.cp2g(samfile.get_reference_name(rd.reference_id),rd.pos+1,C) #self
    
      if covVec[int(gip_GoodReadLeft/100)] > cov_cut \
         or covVec[int(gip_GoodReadLeft/100)-1] > cov_cut \
         or covVec[int(gip_GoodReadLeft/100)+1] > cov_cut:
        continue #NOTE: skip if encouting reads from high coverage region or neighbors
    
      N_goodLeftRead += 1   #number of good read pairs within N_qc
    
      if N_goodLeftRead < statReadMax: #on the same chromosome
        rpinsert.append(rd.template_length)
        readlen.append(rd.infer_query_length())
    
      if debug:
        assert (samfile.get_reference_name(rd.reference_id),rd.pos+1) == lib10x.g2cp(gip_GoodReadLeft,C,acclen_C,chrname), \
        "".join([ "to=", str((samfile.get_reference_name(rd.reference_id),rd.pos+1)), ",gip=", str(gip_GoodReadLeft), ",back=", str(lib10x.g2cp(gip_GoodReadLeft,C,acclen_C,chrname)), "read=", str(rd)])
    
      for alpha in setAlpha: #barfragdata[alpha]=dict(b:{fb1,fb2,...}),fragdata[alpha]=[f1,f2,...],N_f[alpha]=counter 
        if not barfragdata[alpha].has_key(bar): #barcode is new
          barfragdata[alpha][bar]=[N_f[alpha]] #add this molecule to the first of bardcoded frag list
          fragdata[alpha].append([gip_GoodReadLeft,bar,[gip_GoodReadLeft]]) #add a new molecule
          N_f[alpha] += 1;  #move molecule counter next
        else: #barcode exists, 
          if (gip_GoodReadLeft - fragdata[alpha][barfragdata[alpha][bar][-1]][-1][-1]) <= alpha and \
              gip_GoodReadLeft < A[samfile.get_reference_name(rd.reference_id)]: 
            #chcking if within dist and not crossing chr border
            #compare this read gip with last same-barcoded read gip, if within alpha => the same molecule
            #if (gip_GoodReadLeft - fragdata[alpha][barfragdata[alpha][bar][-1]][-1][-1]) <= 0: 
            if (gip_GoodReadLeft - fragdata[alpha][barfragdata[alpha][bar][-1]][-1][-1]) <= 0:
              #DEBUG: stop if added matches tail
              #print >>sys.stderr, "!!! Reversed order is impossible in a sorted bam:"
              #print >>sys.stderr, fragdata[alpha][barfragdata[alpha][bar][-1]], chrname[rd.rname], rd.pos, gip_GoodReadLeft
              continue
            fragdata[alpha][barfragdata[alpha][bar][-1]][-1].append(gip_GoodReadLeft) # append read to this frag's read list
          else: #start an new molecule
            barfragdata[alpha][bar].append(N_f[alpha]) #add new molecule to same barcoded frag list
            fragdata[alpha].append([gip_GoodReadLeft,bar,[gip_GoodReadLeft]]) #add new molecule to fragdata
            N_f[alpha] += 1
        #DEBUG: make sure head is always valid compare to tail
        #cp1 = lib10x.g2cp(fragdata[alpha][N_f[alpha]-1][-1][0],C,acclen_C,chrname)
        #cp2 = lib10x.g2cp(fragdata[alpha][N_f[alpha]-1][-1][-1],C,acclen_C,chrname)
        #assert cp1[0] == cp2[0] and cp1[1] <= cp2[1], "incorrect frag gip tip %d (%s,%d) and end %d (%s,%d) of %s" % (fragdata[alpha][frag][-1][0],cp1[0],cp1[1],fragdata[alpha][frag][-1][-1],cp2[0],cp2[1],str(fragdata[alpha][frag]))
    
  else:  #end of longRangerVersion
  
    print >>sys.stderr, "using 10x longranger inferred molecule..."
    # need to merging in parsing for older longRangerVersion to this one
  
    setAlpha = [0]
    numAlpha = len(setAlpha)
    barfragdata = dict(zip(setAlpha,[ dict() for x in xrange(numAlpha)])) #NOTE: * doesn't create independent dictionaries 
    fragdata = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    N_f =  dict(zip(setAlpha,[0]*numAlpha)) #alpha is alpha, here * does create independent integers
    frag2mol = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    frag2hp = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    frag2ps = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    frag2bar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
  
    for rd in samfile.fetch():
  
      if N_all % 100000 == 0 and N_all != 0: print >>sys.stderr, N_all, "reads processed"
  
      if N_all >= maxRead: break
    
      N_all += 1
  
      if not dict(rd.tags).get('MI'): 
        #throw Molecule unresolved read pairs
        continue
  
      #filter based on sam flags
        # i.e. -f 0x001 -F 0x004 -F 0x008 -F 0x100 -F 0x800 -F 0x400 -F 0x200
        #throw unmapped and low quality read and qc failed read read pairs
      if (not rd.is_paired) or rd.is_unmapped or rd.mate_is_unmapped or rd.is_qcfail or rd.is_secondary or rd.is_supplementary or rd.is_duplicate:
        continue
  
      #filter based on read quality
      if rd.mapq<30:
        continue
  
      try:
         hp = int(dict(rd.tags)['HP'])
         ps = int(dict(rd.tags)['PS'])
      except:
         hp = -1
         ps = -1
    
      #output imBedpe/ffBedpe/rrBedpe

      #include only one ordered pair from two reads
      reference_name, reference_start, next_reference_name, next_reference_start, ordered = lib10x.orderedPair(rd)
      if ( not rd.is_reverse ) and rd.mate_is_reverse and ( not rd.is_proper_pair ) : #IM
        imBedpeFile.write('\t'.join([reference_name, str(reference_start), str(reference_start+100), \
                         next_reference_name, str(next_reference_start), str(next_reference_start+100), str(hp)])+'\n' )
  
      if ( not rd.is_reverse ) and ( not rd.mate_is_reverse ) and ordered : #FF
        ffBedpeFile.write('\t'.join([reference_name, str(reference_start), str(reference_start+100), \
                         next_reference_name, str(next_reference_start), str(next_reference_start+100), str(hp)])+'\n' )
  
      if rd.is_reverse and rd.mate_is_reverse and ordered: #RR
        rrBedpeFile.write('\t'.join([rd.reference_name, str(rd.reference_start), str(rd.reference_start+100), \
                         rd.next_reference_name, str(rd.next_reference_start), str(rd.next_reference_start+100), str(hp)])+'\n' )
  
      N_qc += 1
  
      if rd.pos>rd.mpos: #we only consider the leftmost read of any read pair
        continue
  
      if rd.pos==rd.mpos and rd.is_read2: #if it cann't be distinguished by pos we distinguish by first/second read
        continue #NOTE: this make sure one read pair is only counted once and reduce excessive zeros in spacing
      #NOTE: pysam pos is 0-based !!!!
  
      if rd.pos<0 or rd.pos>=A[samfile.get_reference_name(rd.reference_id)]:
        continue #NOTE: if exceeding the bounds of reference length ignore
  
      bar = dict(rd.tags)['BX'] #BX is the barcode field
      mol = int(dict(rd.tags)['MI']) #MI is the molecule index
  
      gip_GoodReadLeft = lib10x.cp2g(samfile.get_reference_name(rd.reference_id),rd.pos+1,C) #self
    
      try:
        if covVec[int(gip_GoodReadLeft/100)] > cov_cut \
          or covVec[int(gip_GoodReadLeft/100)-1] > cov_cut \
          or covVec[int(gip_GoodReadLeft/100)+1] > cov_cut:
          continue #NOTE: skip if encouting reads from high coverage region or neighbors
      except:
        continue
     
      N_goodLeftRead += 1   #number of good read pairs within N_qc
  
      if N_goodLeftRead < statReadMax: #on the same chromosome
        rpinsert.append(rd.template_length)
        readlen.append(rd.infer_query_length())
  
      # relie on MI
      for alpha in setAlpha: #barfragdata[alpha]=dict(b:{fb1,fb2,...}),fragdata[alpha]=[f1,f2,...],N_f[alpha]=counter
        if not barfragdata[alpha].has_key(bar): #barcode is new
          barfragdata[alpha][bar]=[N_f[alpha]] #add this molecule to the first of bardcoded molecule list
          frag2mol[alpha].append(mol)
          frag2hp[alpha].append(hp)
          frag2ps[alpha].append(ps)
          frag2bar[alpha].append(bar)
          fragdata[alpha].append([gip_GoodReadLeft,bar,[gip_GoodReadLeft]]) #add a new molecule
          N_f[alpha] += 1;  #move molecule counter next
        else: #barcode exists,
          if mol == frag2mol[alpha][barfragdata[alpha][bar][-1]]: #from a registered molecule, expand it
            fragdata[alpha][barfragdata[alpha][bar][-1]][-1].append(gip_GoodReadLeft) #add read to molecule
          else: #from another molecule of the same bar, register it
            frag2mol[alpha].append(mol)
            frag2hp[alpha].append(hp)
            frag2ps[alpha].append(ps)
            frag2bar[alpha].append(bar)
            fragdata[alpha].append([gip_GoodReadLeft,bar,[gip_GoodReadLeft]]) #add a new molecule
            barfragdata[alpha][bar].append(N_f[alpha]) #add this molecule to the list of bardcoded molecule
            N_f[alpha] += 1;  #move molecule counter next
        assert N_f[alpha] == len(fragdata[alpha]) #making sure counting and expanding is consistent
  #      if mol == 1196356 and debug:
          #print(bar)
          #print(barfragdata[alpha][bar])
          #print(frag2mol[alpha])
  #        print("rd.pos+1=", rd.pos+1, " gip=", gip_GoodReadLeft)
          
          #print(fragdata[alpha][barfragdata[alpha][bar][-1]]) #newly added mol
  
  samfile.close()
  imBedpeFile.close()
  ffBedpeFile.close()
  rrBedpeFile.close()
    
  #DEBUG:
  if debug:
    print >>sys.stderr, "N_f[alpha]=", N_f
    print >>sys.stderr, "len(fragdata[alpha])=", dict(zip(setAlpha, [len(fragdata[alpha]) for alpha in setAlpha]))
  #for alpha in setAlpha:
  #  assert N_f[alpha] == len(fragdata[alpha])
  #  for fi in range(0,N_f[alpha]):
  #    cp1 = lib10x.g2cp(fragdata[alpha][fi][-1][0],C,acclen_C,chrname)
  #    cp2 = lib10x.g2cp(fragdata[alpha][fi][-1][-1],C,acclen_C,chrname)
  #    assert cp1[0] == cp2[0] and cp1[1] <= cp2[1], "incorrect frag gip tip %d (%s,%d) and end %d (%s,%d) of %s" % (fragdata[alpha][fi][-1][0],cp1[0],cp1[1],fragdata[alpha][fi][-1][-1],cp2[0],cp2[1],str(fragdata[alpha][fi]))
  
    #well, so right now we are adding in any read instead of read pair 
    #in deriving fragdata and barfragdata
    #so we need to divide by 2 when assessing distribution S(l): fragBarFragRdpCount
    #and also need to remove the effect of neighboring r1 and r2 for D(d): multiRdpFragRdpDistance
  
  #DEBUG
  #for f in barfragdata[50000]['AACACCTCAATCTA-1']:
  #  print >>sys.stderr, fragdata[50000][f] 
  
  #  print >>sys.stderr
    print >>sys.stderr, maxRead, "reads analyzed"
    print >>sys.stderr, "medianInsert=", numpy.median(rpinsert)
    print >>sys.stderr, "madInsert=", lib10x.mad(rpinsert)
    print >>sys.stderr, "medianReadLength=", numpy.median(readlen)
  print >>sys.stderr, N_qc, "of", N_all, "reads passed qc and have BX field"
  print >>sys.stderr, N_goodLeftRead,"read pairs of", N_qc, "reads are good read one with barcodes"
  #print >>sys.stderr, "found", N_b, "total barcodes"
  #print >>sys.stderr
  
  #4.0 summary statistics of barcoded read and molecule
  #brcount=[len(barrpdata[b]["c1"]) for b in barrpdata]
  #print >>sys.stderr, brcount
  #N_fm=copy.deepcopy(N_f)
  #for x in N_fm.keys():  N_fm[x] -= N_bs
  #print >>sys.stderr, "given alpha, found", N_fm, "multi-read molecules for", N_bm, "multi read barcodes estimated"
  #4.0.1 Let's decide which alpha to choose
  #4.0.2 do we want to exclude single appearance barcodes (guess yes)
  #quit()
  
  if True:
    #print >>sys.stderr
    print >>sys.stderr, "performing barcode-molecule analysis..."
    #fragBarBedfile = dict(zip(setAlpha,[re.sub("$",".fragBar.bed",outPrefix)]*numAlpha))
    #for alpha in setAlpha: fragBarBedfile[alpha]=open(".".join([fragBarBedfile[alpha],"a"+str(alpha)]),'w')
    #fragBarTagfile = dict(zip(setAlpha,[re.sub("$",".fragBar.tag",outPrefix)]*numAlpha))
    #for alpha in setAlpha: fragBarTagfile[alpha]=open(".".join([fragBarTagfile[alpha],"a"+str(alpha)]),'w')
  
    multiFragBarBedfile = dict(zip(setAlpha,[re.sub("$",".multiFragBar.bed",outPrefix)]*numAlpha))
    for alpha in setAlpha: multiFragBarBedfile[alpha]=open(".".join([multiFragBarBedfile[alpha],"a"+str(alpha)]),'w')
    multiFragBarTagfile = dict(zip(setAlpha,[re.sub("$",".multiFragBar.tag",outPrefix)]*numAlpha))
    for alpha in setAlpha: multiFragBarTagfile[alpha]=open(".".join([multiFragBarTagfile[alpha],"a"+str(alpha)]),'w')
  
    #print >>sys.stderr, "bardatafile=", bardatafile
    #print >>sys.stderr, "mbardatafile=", mbardatafile
    #print >>sys.stderr, "fragdatafile=", fragdatafile
    #print >>sys.stderr, "multiFragdatafile=", multiFragdatafile
  
    print >>sys.stderr, "chosen alpha=", setAlpha
    #initialize empty dicts
    #fragBar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    multiFragBar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    #fragBarMap = dict(zip(setAlpha,[ dict() for x in xrange(numAlpha)]))
    multiFragBarMap = dict(zip(setAlpha,[ dict() for x in xrange(numAlpha)]))
    #fragCountFragBar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    fragCountMultiFragBar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    #fragLengthFragBar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    fragLengthMultiFragBar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
    
    #for alpha in setAlpha: 
    #  fragBar[alpha] = list(set([ b for b in barfragdata[alpha] ])) #total barcode space, ordered
    #  fragBarMap[alpha] = dict(zip(fragBar[alpha],xrange(1,len(fragBar[alpha])+1))) #barcode to barcode index
    #print >>sys.stderr, "given alpha, found", dict(zip(setAlpha, [ len(fragBar[alpha]) for alpha in setAlpha ])), "total barcodes (fragBar)"
  
    for alpha in setAlpha: 
      multiFragBar[alpha] = list(set([ b for b in barfragdata[alpha] if len(barfragdata[alpha][b])>1 ])) #reduce to barcode having >=2 molecules
      multiFragBarMap[alpha] = dict(zip(multiFragBar[alpha],xrange(0,len(multiFragBar[alpha])))) #0-based index of tags within multiFragBar 
    print >>sys.stderr, "given alpha, found", dict(zip(setAlpha, [ len(multiFragBar[alpha]) for alpha in setAlpha ])), "total barcodes (multiFragBar)"
  
    #for alpha in setAlpha: 
    #  fragCountFragBar[alpha] = [ len(barfragdata[alpha][b]) for b in fragBar[alpha] ]
    #print >>sys.stderr, "given alpha, found", dict(zip(setAlpha, [ sum(fragCountFragBar[alpha]) for alpha in setAlpha ])), "total molecules (fragBar)"
  
    for alpha in setAlpha: 
      fragCountMultiFragBar[alpha] = [ len(barfragdata[alpha][b]) for b in multiFragBar[alpha] ]
    print >>sys.stderr, "given alpha, found", dict(zip(setAlpha, [ sum(fragCountMultiFragBar[alpha]) for alpha in setAlpha ])), "total molecules (multiFragBar)"
  
    #for alpha in setAlpha:
    #  for bar in fragBar[alpha]: 
    #    for frag in barfragdata[alpha][bar]:
    #      cp1 = lib10x.g2cp(fragdata[alpha][frag][-1][0],C,acclen_C,chrname)
    #      cp2 = lib10x.g2cp(fragdata[alpha][frag][-1][-1],C,acclen_C,chrname)
    #      #DEBUG: testing assertion
    #      #assert cp1[0] == cp2[0] and cp1[1] <= cp2[1], "incorrect frag gip tip %d (%s,%d) and end %d (%s,%d) of %s" % (fragdata[alpha][frag][-1][0],cp1[0],cp1[1],fragdata[alpha][frag][-1][-1],cp2[0],cp2[1],str(fragdata[alpha][frag]))
    #      if cp1[0] == cp2[0] and cp1[1] <= cp2[1]:
    #        print >>fragBarBedfile[alpha], "\t".join([str(cp1[0]), str(cp1[1]), str(cp2[1]), str(fragBarMap[alpha][bar]), str(len(fragdata[alpha][frag][-1])), str(cp2[1]-cp1[1]), str(len(barfragdata[alpha][bar])), str(frag2mol[alpha][frag]), str(frag2hp[alpha][frag]), str(frag2ps[alpha][frag]) ])
    ##chr, start, end, barcode_index(0-based), molecule_size(# of reads), molecule_len(base pairs), barcode_duplicates, mi, hp, ps
    #        fragLengthFragBar[alpha].append(cp2[1]-cp1[1])
    #    print >>fragBarTagfile[alpha], bar, len(barfragdata[alpha][bar])
    #    #outputs inferred molecules and its barcode
  
    for alpha in setAlpha:
      for bar in multiFragBar[alpha]: 
        for frag in barfragdata[alpha][bar]:
          cp1 = lib10x.g2cp(fragdata[alpha][frag][-1][0],C,acclen_C,chrname)
          cp2 = lib10x.g2cp(fragdata[alpha][frag][-1][-1],C,acclen_C,chrname)
          #DEBUG: testing assertion
          #assert cp1[0] == cp2[0] and cp1[1] <= cp2[1], "incorrect frag gip tip %d (%s,%d) and end %d (%s,%d) of %s" % (fragdata[alpha][frag][-1][0],cp1[0],cp1[1],fragdata[alpha][frag][-1][-1],cp2[0],cp2[1],str(fragdata[alpha][frag]))
          if cp1[0] == cp2[0] and cp1[1] <= cp2[1]:
            print >>multiFragBarBedfile[alpha], "\t".join([str(cp1[0]), str(cp1[1]), str(cp2[1]), str(multiFragBarMap[alpha][bar]), str(len(fragdata[alpha][frag][-1])), str(cp2[1]-cp1[1]), str(len(barfragdata[alpha][bar])), str(frag2mol[alpha][frag]), str(frag2hp[alpha][frag]), str(frag2ps[alpha][frag]), str(frag2bar[alpha][frag]) ])
            fragLengthMultiFragBar[alpha].append(cp2[1]-cp1[1])
        print >>multiFragBarTagfile[alpha], bar, len(barfragdata[alpha][bar])
        #outputs inferred molecules and its barcode
  
    #meanFragLengthFragBar=dict(zip(setAlpha, [ numpy.mean(fragLengthFragBar[alpha]) for alpha in setAlpha ]))
    #for alpha in setAlpha:
    #  print >>sys.stderr, "mean multi-frag length (fragBar):", meanFragLengthFragBar
  
    meanFragLengthMultiFragBar=dict(zip(setAlpha, [ numpy.mean(fragLengthMultiFragBar[alpha]) for alpha in setAlpha ]))
    if debug:
      for alpha in setAlpha:
        print >>sys.stderr, "mean multi-frag length (multiFragBar):", meanFragLengthMultiFragBar
      print >>sys.stderr, "elapsed time", time.time()-start_time, "seconds"
    
  # 0. what we have
  #    0.0 fragdata  (yes)
  #    0.1 barfragdata (yes)
  #    0.2 fragBar (yes)
  #    0.3 multiFragBar (yes)
  #    0.4 fragCountFragBar (no)
  #    0.5 fragCountMultiFragBar (yes)
  #    0.6 mbrcount (no)
  #    0.7 mbrspace (no)
  #    what we need
  #    0.8 fragBarFragRdpCount [H]
  #    0.9 multiRdpFragRdpDistance [G]
  #    0.10 fragFracMultiFragBar [T]
  fragFracMultiFragBar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
  fragFracFragBar = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
  #fragBarFragRdpCount = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
  multiFragBarFragRdpCount = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
  #fragBarMultiRdpFragRdpDistance = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
  multiFragBarMultiRdpFragRdpDistance = dict(zip(setAlpha,[ [] for x in xrange(numAlpha)]))
  
  #for alpha in setAlpha:  #
  #  fragBarFragRdpCount[alpha] = [ len(fragdata[alpha][f][-1]) for b in fragBar[alpha] for f in barfragdata[alpha][b] ]
  
  for alpha in setAlpha:  #
    multiFragBarFragRdpCount[alpha] = [ len(fragdata[alpha][f][-1]) for b in multiFragBar[alpha] for f in barfragdata[alpha][b] ]
  
  #for alpha in setAlpha: 
  #  fragBarMultiRdpFragRdpDistanceTemp = [ numpy.diff(fragdata[alpha][f][-1]) for b,lf in barfragdata[alpha].iteritems() for f in lf if len(fragdata[alpha][f][-1])>1 ]
  #  fragBarMultiRdpFragRdpDistance[alpha] = [ x for l in fragBarMultiRdpFragRdpDistanceTemp for x in l ] #collapsing list of lists
  
  for alpha in setAlpha: 
    multiFragBarMultiRdpFragRdpDistanceTemp = [ numpy.diff(fragdata[alpha][f][-1]) for b,lf in barfragdata[alpha].iteritems() for f in lf if len(fragdata[alpha][f][-1])>1 ]
    multiFragBarMultiRdpFragRdpDistance[alpha] = [ x for l in multiFragBarMultiRdpFragRdpDistanceTemp for x in l ] #collapsing list of lists
  
  for alpha in setAlpha: 
    multiFragBarTa=sum(fragCountMultiFragBar[alpha])
    fragFracMultiFragBar[alpha] = dict(zip(multiFragBar[alpha],[float(x)/multiFragBarTa for x in fragCountMultiFragBar[alpha]]))
  
  #for alpha in setAlpha: 
  #  fragBarTa=sum(fragCountFragBar[alpha])
  #  fragFracFragBar[alpha] = dict(zip(fragBar[alpha],[float(x)/fragBarTa for x in fragCountFragBar[alpha]]))
  
  # 1. Learn the density functions :
  
  # write TenXdat.alpha.h5
  #for alpha in setAlpha:
  #  TenXdat = pandas.HDFStore(re.sub("$",".fragBar.h5.a%s" % alpha,outPrefix))
  #  #    1.1 H(l)
  #  TenXdat['S'] = pandas.Series(lib10x.empirical_density(numpy.array(fragBarFragRdpCount[alpha])))
  #  #    1.2 G(d)
  #  TenXdat['D'] = pandas.Series(lib10x.empirical_density(numpy.array(fragBarMultiRdpFragRdpDistance[alpha])))
  #  #    1.3 T(t)
  #  TenXdat['T'] = pandas.DataFrame([fragFracFragBar[alpha]])
  #  #    1.4 alpha
  #  TenXdat['alpha'] = pandas.DataFrame([alpha])
  #  #     1.5 N_f
  #  TenXdat['N_f'] = pandas.DataFrame([N_f[alpha]])
  #  #     1.6 N_b
  #  TenXdat['L_f'] = pandas.DataFrame([meanFragLengthFragBar[alpha]])
  #  # Have to be a pandas object to get stored
  #  TenXdat.close()
  #  print >>sys.stderr, "done saving", str(TenXdat), "; elapsed time", time.time()-start_time, "seconds"
  
  for alpha in setAlpha:
    TenXdat = pandas.HDFStore(re.sub("$",".multiFragBar.h5.a%s" % alpha,outPrefix))
    #    1.1 H(l)
    TenXdat['S'] = pandas.Series(lib10x.empirical_density(numpy.array(multiFragBarFragRdpCount[alpha])))
    #    1.2 G(d)
    TenXdat['D'] = pandas.Series(lib10x.empirical_density(numpy.array(multiFragBarMultiRdpFragRdpDistance[alpha])))
    #    1.3 T(t)
    TenXdat['T'] = pandas.DataFrame([fragFracMultiFragBar[alpha]])
    #    1.4 alpha
    TenXdat['alpha'] = pandas.DataFrame([alpha])
    #     1.5 N_f
    TenXdat['N_f'] = pandas.DataFrame([N_f[alpha]])
    #     1.6 N_b
    TenXdat['L_f'] = pandas.DataFrame([meanFragLengthMultiFragBar[alpha]])
    # Have to be a pandas object to get stored
    TenXdat.close()
    if debug:
      print >>sys.stderr, "done saving", str(TenXdat), "; elapsed time", time.time()-start_time, "seconds"
  
  #covChr=covTable.iloc[0,0]
  #covChrIdx=0
  #start_time=time.time()
  #for i in xrange(0,len(covTable)+1):
  #  if debug and i % 10000 == 0:
  #    print >>sys.stderr, "processed", i, "of", len(covTable), "pct", float(i)/len(covTable)
  #    print >>sys.stderr, "remaining time", (time.time()-start_time)/(i+1)*len(covTable)
  #  if i == len(covTable) or covTable.iloc[i,0] != covChr:
  #    covEnd=int(lib10x.cp2g(covTable.iloc[i-1,0],covTable.iloc[i-1,1],C)/100)
  #    if debug:
  #      assert (covEnd-covStart)==((i-1)-covChrIdx)
  #    covVec[covStart:(covEnd+1)] = covTable.iloc[covChrIdx:i,3]
  #    if i == len(covTable): break
  #    covChrIdx=i
  #    covChr=covTable.iloc[i,0]
  #    covStart=int(lib10x.cp2g(covTable.iloc[i,0],covTable.iloc[i,1],C)/100)
  #covEnd=int(lib10x.cp2g(covTable.iloc[len(covTable)-1,0],covTable.iloc[len(covTable)-1,1],C)/100)
  #covVec[covStart:(covEnd+1)] = covTable.iloc[:,3]
  #for i in xrange(0,len(covTable)):
  #  covVec[int(lib10x.cp2g(covTable.iloc[i,0],covTable.iloc[i,1],C)/100)] = covTable.iloc[i,3]

#with open(barrpdatafile, 'w') as outfile:  json.dump(barrpdata, outfile)
#with open(refdatafile, 'w') as outfile:  json.dump(refdata, outfile)

#can remove rpinsert and readlen
#80G mem OK up to here

#4.0 right now the only way of exchanging data between R and Python is simple text output
#none of these data files are used, let's skip creating them later?
#if False: #this is currently disabled
#  print >>sys.stderr
#  print >>sys.stderr, "performing barcode-read analysis..."
#  
#  brcountdat = open(re.sub("$",".brcount.dat",bamfile),"w")
#  for x in brcount:  print >>brcountdat, x
#  N_bm=sum([x>1 for x in brcount])
#  del brcount; print >>sys.stderr, "recycle brcount: ", gc.collect()
#  
#  brspace = [numpy.diff(barrpdata[b]["p_GoodReadLeft"])[numpy.diff(barrpdata[b]["c1"])==0] for b in barrpdata] #list
#  brspacedat = open(re.sub("$",".brspace.dat",bamfile),"w")
#  for x in [ x for s in brspace for x in s ]:  print >>brspacedat, x
#  del brspace; print >>sys.stderr, "recycle brspace: ", gc.collect()
#  
#  brinsert = [ abs(numpy.array([x - y for x, y in zip(barrpdata[b]["p_GoodReadLeft"],barrpdata[b]["p2"])])[numpy.equal(barrpdata[b]["c1"],barrpdata[b]["c2"])]) for b in barrpdata]
#  brinsertdat = open(re.sub("$",".brinsert.dat",bamfile),"w")
#  for x in [ x for s in brinsert for x in s ]:  print >>brinsertdat, x
#  del brinsert; print >>sys.stderr, "recycle brinsert: ", gc.collect()
#  
#  barcode = [ b for b in barrpdata.keys() ] #br_list
#  print >>sys.stderr, "barcode-read analysis found", len(barcode), "total barcodes"
#  del barcode; print >>sys.stderr, "recycle barcode: ", gc.collect()
#  
#  mbarcode = [ b for b in barrpdata if len(barrpdata[b]["c1"])>1 ] #mbr_list
#  print >>sys.stderr, "barcode-read analysis found", len(mbarcode), "multi-read barcodes"
#  
#  mbrcount = [ len(barrpdata[b]["c1"]) for b in mbarcode ]
#  mbrcountdat = open(re.sub("$",".mbrcount.dat",bamfile),"w")
#  for x in mbrcount:  print >>mbrcountdat, x
#  del mbrcount; print >>sys.stderr, "recycle mbrcount: ", gc.collect()
#  
#  mbrspace = [ numpy.diff(barrpdata[b]["p_GoodReadLeft"])[numpy.diff(barrpdata[b]["c1"])==0] for b in mbarcode ]
#  mbrspacedat = open(re.sub("$",".mbrspace.dat",bamfile),"w")
#  for x in [ x for s in mbrspace for x in s ]:  print >>mbrspacedat, x
#  del mbrspace; print >>sys.stderr, "recycle mbrspace: ", gc.collect()
#  
#  mbrinsert = [ abs(numpy.array([x - y for x, y in zip(barrpdata[b]["p_GoodReadLeft"],barrpdata[b]["p2"])])[numpy.equal(barrpdata[b]["c1"],barrpdata[b]["c2"])]) for b in mbarcode]
#  mbrinsertdat = open(re.sub("$",".mbrinsert.dat",bamfile),"w")
#  for x in [ x for s in mbrinsert for x in s ]:  print >>mbrinsertdat, x
#  del mbrinsert; print >>sys.stderr, "recycle mbrinsert: ", gc.collect()
#
#
#  N_b=len(barrpdata)
#  N_bs=N_b-N_bm
#  print >>sys.stderr, "barcode-read analysis:", N_bm, "of", N_b, "barcodes appeared more than one times"
#  print >>sys.stderr, "barcode-read analysis numbers may seem contradictory to molecule analysis numbers in small region"
#  print >>sys.stderr, "these numbers are expected to be close to each other in global"
#  print >>sys.stderr

#barrpdatafile = re.sub("$",".barrp.json",bamfile)
#refdatafile = re.sub("$",".ref.json",bamfile)
#hdfdatafile = pandas.HDFStore(re.sub("bam$","data.hd5",bamfile))


def run():
  parser = argparse.ArgumentParser()
  parser.add_argument("bam",help="sorted 10x bamfile")
  parser.add_argument("-c","--covFile", default="none",help="100bp window coverage file")
  parser.add_argument("-l","--longRangerVersion", type=float, default=1.2,
    help="v>1.2 will use 10X's MI tag to infer molecules if specified and available")
  parser.add_argument("-f","--setAlpha",default="30000,50000,70000",
    help="alpha i.e. max molecule extension basepair, override by 0 in case of using MI tag")
  parser.add_argument("-b","--maxRead",type=int,default=numpy.iinfo(numpy.int64).max,help="max number of read to process")
  parser.add_argument("-o","--outPrefix",default="none",help="outputPrefix")
  parser.add_argument("-a","--debug",default=False,action="store_true",help="in debug mode")
  try:
    args = parser.parse_args()
  except SystemExit:
    args = None
  main(args)
  print >>sys.stderr, "PARSE10XDONE; elapsed time", time.time()-start_time, "seconds"
  
if __name__ == '__main__':
  run()
