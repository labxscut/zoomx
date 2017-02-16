#!/usr/bin/env python
#using two spaces as indent
#this script assuming bam is sorted by coordinate
import argparse
import os
import sys
import pandas
import pybedtools
import pysam
import lib10x
import time
import re
import numpy
import collections
from multiprocessing import Pool

"""
Define parameters
=================
"""
parser = argparse.ArgumentParser()
parser.add_argument("bam",help="sorted 10x bamfile")
parser.add_argument("-e","--empirical",help=".h5 file contains empricial distribution")
parser.add_argument("-r","--region",help=".bed/.bedpe file defines scanning regions")
parser.add_argument("-f","--regformat",help=".bed/.bedpe")
parser.add_argument("-d","--delta",help="size of scan trunk")
parser.add_argument('-n', '--nprocs', dest='nprocs', default=1, type=int, help="split into multiple processes (default: %(default)s )")
parser.add_argument("-x","--insertExtend",type=int,default=200, help="bp extension to cover break points")
parser.add_argument("-w","--windowStep",type=int,default=10, help="bp stepsize to create region windows")

try:
	args = parser.parse_args()
	bamfile = args.bam
	exmpirical = args.empirical
	regionFile = args.region
	regformat = args.regformat
	insertExtend = args.insertExtend
	windowStep = args.windowStep
	delta = args.delta
	nprocs = args.nprocs
except SystemExit:
	bamfile = os.environ['HOME']+"/ws/bcells/P04979.wg10.bam"
	empiricalfile = bamfile + ".a50000.h5"
	regionFile = bamfile + ".bedpe.filt.new.head"
	regformat = "bedpe"
	insertExtend = 200
	windowStep = 10
	delta = 50000
	nprocs = 1

start_time=time.time()
#refdata = lib10x.load_refdata(bamfile)

"""
2. trunk u and v calculation
============================
"""

"""
	2.0 Split scan regions into delta spaced intervals as BedTools Objects
	------------------------------------------------------------------------
			2.0.1 Create a sorted region list: allRegion (1 ... m)

"""
#NOTE:
#	if input is bed, create sliding windows cover pairwise combination
#	allRegion: containing all relavant regions; [ 
#	allLogLR: containing all index of region pairs to be scanned
allRegion = [] # pandas.dataframe
allLogLR = []	# list of two tuples
regionIndex = 0
logLRIndex = 0
if regformat == 'bed':
	allRegion = pybedtools.BedTool(regionFile)
	allRegion = allRegion.window_maker(b=allRegion,w=2*delta) # 351 Regions of Chr22
	allRegion = allRegion.to_dataframe()
	allLogLR = [ (x,y) for x in range(len(allRegion)-1) for y in range(x+1,len(allRegion)) ]
elif regformat == 'bedpe':
	# RATIONAL:
	#	 regions always start from its left limit: [u-d, u, u+d]
	#	 if a read is forward strand (1), break is next to its right upto insertExtend distance
	#	 if a read is backward strand (0), break is next to its left upto insertExtend distance
	#	 we take a bedpe pair
	#	 determine the break point is to right or left
	#	 determine the possible range of break point (insertExtend) toward right or left
	#	 determine the window step size
	#	 generate and transform the windows cover the possible breakpoint range into standard format 
	#	 generate pairwise combination of windows indexes
	#		append windows to allRegions
	#	 append combinations to allLogLR	
	labels = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name','score', 'strand1', 'strand2'] \
							+ ['ft_'+bamfile,'bg_'+bamfile] \
							+ ['note']
	print >>sys.stderr, "labels=", labels
	bedpeInput = pybedtools.BedTool(regionFile).to_dataframe(comment='#',names=labels,low_memory=False,sep="\t")
	bed1 = bedpeInput.iloc[:,[0,1,2,8]] 
	bed2 = bedpeInput.iloc[:,[3,4,5,9]]
	for idx in range(len(bed1)):
		if bed1.iloc[idx,3] == 1:
			br1 = numpy.arange(bed1.iloc[idx,2], bed1.iloc[idx,2] + insertExtend, windowStep)
		else:
			br1 = numpy.arange(bed1.iloc[idx,1] - insertExtend, bed1.iloc[idx,1], windowStep)
		if bed2.iloc[idx,3] == 1: 
			br2 = numpy.arange(bed2.iloc[idx,2], bed2.iloc[idx,2] + insertExtend, windowStep)
		else:
			br2 = numpy.arange(bed2.iloc[idx,1] - insertExtend, bed2.iloc[idx,1], windowStep)
		br1Region = pandas.DataFrame({"chrom":bed1.iloc[idx,0], "start":br1-delta, "end":br1+delta}, columns=["chrom","start","end"])
		br2Region = pandas.DataFrame({"chrom":bed2.iloc[idx,0], "start":br2-delta, "end":br2+delta}, columns=["chrom","start","end"])
		allRegion.append(br1Region)
		allRegion.append(br2Region)
		regionIndex += len(br1Region)+len(br2Region)
		allLogLR.append( [(logLRIndex+x,logLRIndex+len(br1Region)+y) for x in range(len(br1Region)) for y in range(len(br2Region))] )
		logLRIndex += len(br1Region) * len(br2Region)
		""" DEBUG:
		print >>sys.stderr, "br1Region=", len(br1Region), ",br2Region=", len(br2Region)
		print >>sys.stderr, br1Region, br2Region
		print >>sys.stderr, allLogLR
		print >>sys.stderr, "regionIndex=", regionIndex, ",logLRIndex=", logLRIndex
		quit()
		"""
	allRegion = pandas.concat(allRegion)
else:
	print >>sys.stderr, "unacceptable reg file format"

Empirical_H = pandas.read_hdf(empiricalfile,'H')
Empirical_G = pandas.read_hdf(empiricalfile,'G')
Empirical_T = pandas.read_hdf(empiricalfile,'T')
Empirical_alpha = pandas.read_hdf(empiricalfile,'alpha')

"""
	2.1 Parallely carry out the computation
	-------------------------------------------------------------------

	2.1.1 Parallelize the pairwise computation by pool and map pairs: [(1,1):, ..., (m-1,m):]

		Assumptions:
			For each	pair of region (u,v) onstructed, putative break point is assumed to be the middle point of u and v.

		Will use a pair indexed dict: allLogLR [(1,1):, ..., (m-1,m):]

"""

#NOTE:
# for dataframe loc is identifying by label
#                iloc is identifying by order
#DEBUG:
#	print >>sys.stderr, allRegion
#	print >>sys.stderr, allLogLR
#	print >>sys.stderr, allRegion['chrom'].iloc[0]
#	print >>sys.stderr, allRegion.iloc[0].loc['chrom']

allLoad = [ [ lib10x.load_quad(bamfile, (allRegion['chrom'].iloc[i], allRegion['start'].iloc[i], int(0.5*(allRegion['start'].iloc[i]+allRegion['end'].iloc[i]))-1), Empirical_T.keys(), i), \
						  lib10x.load_quad(bamfile, (allRegion['chrom'].iloc[i], int(0.5*(allRegion['start'].iloc[i]+allRegion['end'].iloc[i]))-1, allRegion['end'].iloc[i]), Empirical_T.keys(), i) ] \
						 for i in range(len(allRegion)) ]

def process_pair(it):
	i, alt, mode, pi = it
	print >>sys.stderr, "processing", i, " ", alt, " ", mode, " ", pi 
	#i, alt, mode, pi = (i,'C',4,1.0)
	uidx, vidx = allLogLR[i]
	urp = allLoad[uidx]
	vrp = allLoad[vidx]
	u=int(0.5*(allRegion['start'][uidx]+allRegion['end'][uidx]))
	v=int(0.5*(allRegion['start'][vidx]+allRegion['end'][vidx]))
	return lib10x.compute_uv(u,v,urp,vrp,Empirical_H,Empirical_G,Empirical_T,alt,mode,pi)

print >>sys.stderr, "total", len(allRegion), "regions"
print >>sys.stderr, "allRegion", allRegion
print >>sys.stderr, "total", len(allLogLR), "comparisons"
print >>sys.stderr, "parallilized to", nprocs, "processes"

"""

NOTE: Testing compute_lrt

LogLR1=[]; LogLR2=[]

if nprocs<2:
	for i in range(len(allLogLR)):
		LogLR1.append(process_pair((i,'C',4,1.0)))
else:
	pool = Pool(processes=int(nprocs))
	it = zip(range(len(allLogLR)),'C'*len(allLogLR),[4]*len(allLogLR),[1.0]*len(allLogLR))
	print >>sys.stderr,it
	LogLR1 = pool.map(process_pair,it)

if nprocs<2:
	for i in range(len(allLogLR)):
		LogLR2.append(process_pair((i,'X',4,1.0)))
else:
	pool = Pool(processes=int(nprocs))
	it = zip(range(len(allLogLR)),'X'*len(allLogLR),[4]*len(allLogLR),[1.0]*len(allLogLR))

	LogLR2 = pool.map(process_pair,it,1)

print(LogLR1)
print(LogLR2)

quit()
"""

"""
LogLR1=[]; LogLR2=[]

if nprocs<2:
	for i in range(len(allLogLR)):
		LogLR1.append(process_pair((i,'C',1,1.0)))
else:
	pool = Pool(processes=int(nprocs))
	it = zip(range(len(allLogLR)),'C'*len(allLogLR),[1]*len(allLogLR),[1.0]*len(allLogLR))
	print >>sys.stderr,it
	LogLR1 = pool.map(process_pair,it)

if nprocs<2:
	for i in range(len(allLogLR)):
		LogLR2.append(process_pair((i,'X',1,1.0)))
else:
	pool = Pool(processes=int(nprocs))
	it = zip(range(len(allLogLR)),'X'*len(allLogLR),[1]*len(allLogLR),[1.0]*len(allLogLR))

	LogLR2 = pool.map(process_pair,it,1)

print(LogLR1)
print(LogLR2)

quit()
"""

"""
LogLR1=[]; LogLR2=[]

if nprocs<2:
	for i in range(len(allLogLR)):
		LogLR1.append(process_pair((i,'C',2,1.0)))
else:
	pool = Pool(processes=int(nprocs))
	it = zip(range(len(allLogLR)),'C'*len(allLogLR),[2]*len(allLogLR),[1.0]*len(allLogLR))
	print >>sys.stderr,it
	LogLR1 = pool.map(process_pair,it)

if nprocs<2:
	for i in range(len(allLogLR)):
		LogLR2.append(process_pair((i,'X',2,1.0)))
else:
	pool = Pool(processes=int(nprocs))
	it = zip(range(len(allLogLR)),'X'*len(allLogLR),[2]*len(allLogLR),[1.0]*len(allLogLR))

	LogLR2 = pool.map(process_pair,it,1)

print(LogLR1)
print(LogLR2)
"""

LogLR1=[]; LogLR2=[]

if nprocs<2:
	for i in range(len(allLogLR)):
		LogLR1.append(process_pair((i,'C',7,1.0)))
else:
	pool = Pool(processes=int(nprocs))
	it = zip(range(len(allLogLR)),'C'*len(allLogLR),[7]*len(allLogLR),[1.0]*len(allLogLR))
	print >>sys.stderr,it
	LogLR1 = pool.map(process_pair,it)

if nprocs<2:
	for i in range(len(allLogLR)):
		LogLR2.append(process_pair((i,'X',7,1.0)))
else:
	pool = Pool(processes=int(nprocs))
	it = zip(range(len(allLogLR)),'X'*len(allLogLR),[7]*len(allLogLR),[1.0]*len(allLogLR))

	LogLR2 = pool.map(process_pair,it,1)

print(LogLR1)
print(LogLR2)

print >>sys.stderr, "can I reach here"


"""
3. make junctions calls
===========================
"""

"""
	3.1 sort by likelihood score
	----------------------------
"""


"""
	2.3 save necessary information of each junction
	---------------------------------------------------------------
"""

allLogLR=dict(zip(allLogLR,zip(LogLR1,LogLR2)))


"""
	3.2 resolve overlapping call structure
	----------------------------------------
"""

"""
	3.3 resolve call heterogeneity
	--------------------------------
"""


"""
	4 simulation
	===============
"""

"""
	4.1 simulate deletion
	------------------------------
		123456789 -> 1245689 (3,7 deleted)
"""

"""
	4.2 simulate duplication
	-------------------------------
		123456789 -> 12334567789 (3,7 duplicated)
"""

"""
	4.3 simulate transposition
	-----------------------------------
		123456789 -> 127456389 (3,7 switched)
"""

"""
	4.4 inversion / insertion
	-----------------------------
"""

"""
	4.5 overall simulation procedure
	---------------------------------
		4.5.1 set a sequence span of gin positions
		4.5.2 create a set of events and work out the altered mapping positions
		4.5.3 calculate the number of fragments and reads given pi, then for each sequence
		4.5.3 generate a fragment location uniform
		4.5.4 generate a fragment barcode by T
		4.5.5 generate the number of reads on fragment by H
		4.5.6 generate the spacing between reads by G
		4.5.7 add the generated m reads to gin ordered dictionary sim_data i.e. { g1:{b1:n1,b2:n2},...,gm:{b1:n1,b2:n2} }
		4.5.8 merge the sv and non-sv sim_data
		4.5.9 load_quad_sim i.e. loads the simulated data as if from bam, fill in necessary bits
		4.5.10 run the caller
		4.5.11 compare calls and events
"""
