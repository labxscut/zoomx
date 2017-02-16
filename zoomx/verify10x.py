#!/usr/bin/env python
import sys,os,pysam
import argparse
import re
import shelve
import numpy
import json
from collections import Counter

#TODO: this is to replace sorted
def radixsort(random_list):
    len_random_list = len(random_list)
    modulus = 10
    div = 1
    while True:
        # empty array, [[] for i in range(10)]
        new_list = [[], [], [], [], [], [], [], [], [], []]
        for value in random_list:
            least_digit = value % modulus
            least_digit /= div
            new_list[least_digit].append(value)
        modulus = modulus * 10
        div = div * 10
        
        if len(new_list[0]) == len_random_list:
            return new_list[0]
        
        random_list = []
        rd_list_append = random_list.append
        for x in new_list:
            for y in x:
                rd_list_append(y)

parser = argparse.ArgumentParser()
parser.add_argument("bam",help="name sorted 10x bamfile")
args = parser.parse_args()

bamfile = args.bam 
pydatafile = re.sub("bam$","shelve",bamfile)
samfile = pysam.AlignmentFile(bamfile, "rb")

rn = None; rp_same_bar = []; rp_aln_num=[]; rp_all_map=[]; rp_same_frag = []; rp_same_pos = []; rp_one_pos = []; rp_consistent = []; rp_insert = []; rp_single = 0
rpmax=numpy.inf 
#rpmax=10000
for rd in samfile.fetch(until_eof=True): #this will miss the last rp!
    if rd.is_unmapped or rd.is_qcfail or rd.mapping_quality<30 or rd.mate_is_unmapped or rd.is_secondary:
        #throw unmapped and low quality read and qc failed read read pairs
        continue
    if not dict(rd.tags).get('BX'):
        #throw Barcode unresolved read pairs
        continue
    if not rn: #getting a new read pair (rp)
        rn = rd.qname; #print "rn1:",rn,rd.qname
        rb = [ dict(rd.tags)['BX'].split('-')[-1] ]
        pos = [ rd.pos ]
        mpos = [ rd.mpos ]
        ump = [ rd.is_unmapped ]
        chrs = [ rd.rname ]
        mchrs = [ rd.mrnm ]
    elif rn == rd.qname: #within the same rp
        rb.append( dict(rd.tags)['BX'].split('-')[-1] )
        pos.append( rd.pos )
        ump.append( rd.is_unmapped )
        chrs.append( rd.rname )
        mpos.append( rd.mpos )
        mchrs.append( rd.mrnm )
        #print "rn2:",rn,rd.qname
        #print pos, chrs, mpos, mchrs
        #quit()
    else: #wrapping up an old rp and starting a new one
        #print "rn3:",rn,rd.qname
        rp_aln_num.append(len(rb))
        same_bar = len(set(rb))==1
        all_map = not any(ump)
        same_frag = False
        same_chr = len(set(chrs))==1 
        one_pos = len(pos)==1 and all_map
        same_pos = len(pos)>1 and max(numpy.diff(sorted(pos)))<86 and all_map
        consistent = set(mpos) == set(pos) and set(mchrs) == set(chrs) or len(pos) == 1
        if same_bar and all_map and same_chr:
          same_frag = len(pos)==1 or max(numpy.diff(sorted(pos)))<30000
          #NOTE: if pos has nan then the rp_sam_frag with contain False
          #TODO: if speed is needed here, this sorted function can be implemented in O(n) in C
        if not consistent:
          print "rp info:", "rn=", rn, "rb=",rb, "ump=",ump,"pos=",pos, "chrs=",chrs, 'mpos=', mpos, 'mchrs', mchrs
        else:
          rp_insert.append(abs(mpos[0]-pos[0]))
          if(len(pos)==1):
               rp_single=rp_single+1
        #  quit()
        rp_same_bar.append(same_bar)
        rp_all_map.append(all_map)
        rp_same_frag.append(same_frag)
        rp_same_pos.append(same_pos)
        rp_one_pos.append(one_pos)
        rp_consistent.append(consistent)
        #new read pair
        rn = rd.qname
        rb = [ dict(rd.tags)['BX'].split('-')[-1] ]#BX is the barcode field
        pos = [ rd.pos ]
        ump = [ rd.is_unmapped ]
        chrs = [ rd.rname ] 
        mpos = [ rd.mpos ]
        mchrs = [ rd.mrnm ]
        if len(rp_same_bar)>rpmax: #read enough rp
            break

print >>sys.stdout, "total rp:", len(rp_same_bar)
print >>sys.stdout, "rp_same_bar:", sum(rp_same_bar), sum(rp_same_bar)/float(len(rp_same_bar))
print >>sys.stdout, "rp_same_pos:", sum(rp_same_pos), sum(rp_same_pos)/float(len(rp_same_bar))
print >>sys.stdout, "rp_one_pos:", sum(rp_one_pos), sum(rp_one_pos)/float(len(rp_same_bar))
print >>sys.stdout, "rp_all_map:", sum(rp_all_map), sum(rp_all_map)/float(len(rp_same_bar))
print >>sys.stdout, "rp_same_frag:", sum(rp_same_frag), sum(rp_same_frag)/float(len(rp_same_bar))
print >>sys.stdout, "rp_consistent:", sum(rp_consistent), sum(rp_consistent)/float(len(rp_same_bar))
print >>sys.stdout, "rp_aln_num:", Counter(rp_aln_num), 
print >>sys.stdout, "rp_insert:", numpy.mean(rp_insert), "+/-", numpy.std(rp_insert)
print >>sys.stdout, "rp_insert:", numpy.median(rp_insert)

#my_shelf = shelve.open(pydatafile,'n') # 'n' for new
#
#for key in ['rp_same_bar','rp_aln_num','rp_same_frag','rp_all_map']:
#    try:
#        my_shelf[key] = globals()[key]
#    except TypeError:
#        #
#        # __builtins__, my_shelf, and imported modules can not be shelved.
#        #
#        print('ERROR shelving: {0}'.format(key))
#my_shelf.close()
