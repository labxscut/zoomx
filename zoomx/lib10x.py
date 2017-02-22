import bisect
import pysam
import pybedtools
import numpy
import collections
import sys
import os
import pandas
import scipy
import string
import tempfile
from multiprocessing import Pool
from scipy.stats import poisson
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix

# Copyright (C) 2008 Simone Leo - CRS4. All Rights Reserved.
# 
# Permission to use, copy, modify, and distribute this software and its
# documentation for educational, research, and not-for-profit purposes, without
# fee and without a signed licensing agreement, is hereby granted, provided
# that the above copyright notice, this paragraph and the following two
# paragraphs appear in all copies, modifications, and distributions. Contact
# CRS4, Parco Scientifico e Tecnologico, Edificio 1, 09010 PULA (CA - Italy),
# +3907092501 for commercial licensing opportunities.
# 
# IN NO EVENT SHALL CRS4 BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
# INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
# THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF CRS4 HAS BEEN ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# CRS4 SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
# HEREUNDER IS PROVIDED "AS IS". CRS4 HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
# SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

"""
Functions for computing the False Discovery Rate (FDR) of a multiple test
procedure. The FDR is the expected proportion of falsely rejected hypotheses.

@see:
  Yoav Benjamini and Yosef Hochberg, I{Controlling the False Discovery Rate:
  A Practical and Powerful Approach to Multiple Testing}. Journal of the Royal
  Statistical Society, Series B (Methodological), Vol. 57, No. 1 (1995), pp.
  289-300.

  Yekutieli, D. and Benjamini, Y., I{Resampling-based false discovery rate
  controlling multiple test procedures for correlated test statistics}. J. of
  Statistical Planning and Inference, 82, pp. 171-196, 1999.

  John D. Storey (2002), I{A direct approach to false discovery rates}. Journal
  of the Royal Statistical Society: Series B (Statistical Methodology) 64 (3),
  479-498.
"""

__docformat__ = 'epytext en'

import operator
import random

def ThreadIdentifySupportRDP(info):
  #sys.stderr.write("in test thread\n")
  bi=info[0]
  numCall=info[1]
  df=info[2]
  tmpBedpeFile=info[3]
  debug=info[4]
  sys.stderr.write("START: procssing "+str(bi+1)+" of "+str(numCall)+" events\n")
  singleBedpeFile=getTemporaryFile('bedpe_')
  tmpOverlapFile=getTemporaryFile('bedpe_')
  df.to_csv(singleBedpeFile,header=False, sep="\t", index=False)
  pair2pair_cmd="bedtools pairtopair -a %s -b %s | cut -f15 >%s" % (singleBedpeFile, tmpBedpeFile, tmpOverlapFile)
  if debug:
    sys.stderr.write("singleBedpeFile='"+str(singleBedpeFile)+"'\n")
    sys.stderr.write("pair2pair_cmd="+str(pair2pair_cmd)+"\n")
  if run_cmd(pair2pair_cmd) == False:
    sys.stderr.write("Error in running pair2pair overlapping..."); quit()
  overlapCalls=pandas.read_table(tmpOverlapFile, names=["bedpe2"], dtype={'bedpe2':numpy.int32})
  sys.stderr.write("END: procssing "+str(bi+1)+" of "+str(numCall)+" events\n")
  support=list(set(list(overlapCalls['bedpe2']+1)))
  return support

def ParallelIdentifySupportRDP(callDF, rdpBedpeFile, callLabels, callTypes, numprocs, debug):
  #sys.stderr.write("in test rdp\n")
  rdpBedpe=pandas.read_table(rdpBedpeFile, names=callLabels, dtype=dict(zip(callLabels,callTypes)))
  sel=rdpBedpe['chr2']=="="
  rdpBedpe.loc[sel,('chr2')]=rdpBedpe.loc[sel,('chr1')]
  rdpBedpe['index']=xrange(0,rdpBedpe.shape[0]) #0-based index
  tmpBedpeFile=getTemporaryFile('bedpe_')
  #if debug: sys.stderr.write("tmpBedpeFile='"+str(tmpBedpeFile)+"'\n")
  rdpBedpe.to_csv(tmpBedpeFile,header=False, sep="\t", index=False)
  numCall=len(callDF)
  info=[ [i, numCall, callDF.loc[i:i,:], tmpBedpeFile, debug] for i in range(0,numCall) ] #0-based call index
  support=[ [] for i in range(0,numCall) ] #1-based merged call index
  pool=Pool(processes=int(numprocs))
  support = pool.map( ThreadIdentifySupportRDP, info)
  return support

def empirical_density(v):
  #values, base = np.histogram(nearest, bins=20, density=1)
  d, b = numpy.histogram(v, bins=range(0,int(numpy.ceil(v.max()))+1), density=1)
  z = d != 0
  for i in range(0,len(z)):
    if z[i] and i>0:
      j=0
      while (i-j-1>=0 and not z[i-j-1]):
        j+=1
      if i-j != 0: #not splitting in for first zeros
        d[(i-j):(i+1)]=d[i]/(j+1) #fill backward
  return d
#v=np.random.exponential(100,100)
#empirical_density(v)

def orderedPair(rd):
  reference_name = rd.reference_name
  reference_start = rd.reference_start
  next_reference_name = rd.next_reference_name
  next_reference_start = rd.next_reference_start
  ordered = True
  if rd.reference_name > rd.next_reference_name:
    reference_name = rd.next_reference_name
    reference_start = rd.next_reference_start
    next_reference_name = rd.reference_name
    next_reference_start = rd.reference_start
    ordered = False
  elif rd.reference_name == rd.next_reference_name and rd.reference_start>rd.next_reference_start:
    reference_start = rd.next_reference_start
    next_reference_start = rd.reference_start
    ordered = False
  return [reference_name, reference_start, next_reference_name, next_reference_start, ordered]

def mad(a, axis=None):
  """
  Compute *Median Absolute Deviation* of an array along given axis.
  """
  med = numpy.median(a, axis=axis)                # Median along given axis
  if axis is None:
    umed = med                              # med is a scalar
  else:
    umed = numpy.expand_dims(med, axis)         # Bring back the vanished axis
  mad = numpy.median(numpy.absolute(a - umed), axis=axis)  # MAD along given axis
  return mad

def randascii(n):
  n=min(max(n,1),32)  #bound it between 1 to 32
  return ''.join([random.choice(string.ascii_letters) for i in range(n)])

def getTemporaryFile(prefix): #this provides a unique temporary file name only
  #add suffix here to avoid gzip error: already has _z suffix -- unchanged
  return os.path.join(tempfile.tempdir,prefix+next(tempfile._get_candidate_names())+"tmp")

def setTemporaryDir(tempdir):
  tempfile.tempdir = tempdir
  try:
    os.mkdir(tempfile.tempdir)
  except Exception:
    sys.stderr.write("tempdir="+str(tempfile.tempdir)+" error!, either it exists or cann't be created! please change it\n")
    quit()
  return True

def rmTemporaryDir(tempdir):
  try:
    #shutil.rmtree(tempdir)
    os.system("rm -rf {}".format(tempdir))
  except Exception:
    sys.stderr.write("tempdir="+str(tempfile.tempdir)+" error!, it cann't be removed! please remove it mannually\n")
    return False 
  return True

def run_cmd(cmd):
    try:
        os.system(cmd)
    except OSError:
        return False
    return True

#def bedpe_pair2pair(

def bh_rejected(pv, threshold):
    """
    Return the list of rejected p-values from C{pv} at FDR level C{threshold}.

    The downside of this approach is that the FDR level must be chosen in
    advance - the L{bh_qvalues} function yields a p-value-like output instead.

    @type pv: list
    @param pv: p-values from a multiple statistical test
    @type threshold: float
    @param threshold: the level at which FDR rate should be controlled

    @rtype: list
    @return: p-values of rejected null hypotheses
    """
    if threshold < 0 or threshold > 1:
        raise ValueError("the threshold must be between 0 and 1")
    if not pv:
        return []
    pv = sorted(pv)
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    m = len(pv)
    for i in xrange(m-1, -1, -1):
        if pv[i] <= float(i+1)*threshold/m:
            return pv[:i+1]
    return []

def bh_qvalues(pv):
    """
    Return Benjamini-Hochberg FDR q-values corresponding to p-values C{pv}.

    This function implements an algorithm equivalent to L{bh_rejected} but
    yields a list of 'adjusted p-values', allowing for rejection decisions
    based on any given threshold.

    @type pv: list
    @param pv: p-values from a multiple statistical test

    @rtype: list
    @return: adjusted p-values to be compared directly with the desired FDR
      level
    """
    if not pv:
        return []
    m = len(pv)
    args, pv = zip(*sorted(enumerate(pv), None, operator.itemgetter(1)))
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    qvalues = m * [0]
    mincoeff = pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in xrange(m-2, -1, -1):
        coeff = m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    return qvalues


def get_pi0(pv, lambdas):
    """
    Compute Storey's C{pi0} from p-values C{pv} and C{lambda}.

    this function is equivalent to::
    
        m = len(pv)
        return [sum(p >= l for p in pv)/((1.0-l) * m) for l in lambdas]
        
    but the above is C{O(m*n)}, while it needs only be C{O(m+n)
    (n = len(lambdas))}

    @type pv: list
    @param pv: B{SORTED} p-values vector
    @type lambdas: list
    @param lambdas: B{SORTED} lambda values vector

    @rtype: list
    @return: estimated proportion of null hypotheses C{pi0} for each lambda
    """
    m = len(pv)
    i = m - 1
    pi0 = []
    for l in reversed(lambdas):
        while i >= 0 and pv[i] >= l:
            i -= 1
        pi0.append((m-i-1)/((1.0-l)*m))
    pi0.reverse()
    return pi0


def storey_qvalues(pv, l=None):
    """
    Return Storey FDR q-values corresponding to p-values C{pv}.

    The main difference between B-H's and Storey's q-values is that the latter
    are weighted by the estimated proportion C{pi0} of true null hypotheses.

    @type pv: list
    @param pv: p-values from a multiple statistical test
    @type l: float
    @param l: lambda value for C{pi0} (fraction of null p-values) estimation

    @rtype: list
    @return: storey q-values corresponding to C{pv}
    """
    if not pv:
        return []
    m = len(pv)
    args, pv = zip(*sorted(enumerate(pv), None, operator.itemgetter(1)))
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")

    if l is None:
        # optimal lambda/pi0 estimation
        lambdas = [i/100.0 for i in xrange(0, 91, 5)]
        n = len(lambdas)
        pi0 = get_pi0(pv, lambdas)
        min_pi0 = min(pi0)
        mse = [0] * n        
        for i in xrange(1, 101):
            # compute bootstrap sample with replacement
            pv_boot = [pv[int(random.random()*m)] for j in xrange(m)]
            pi0_boot = get_pi0(sorted(pv_boot), lambdas)
            for j in xrange(n):
                mse[j] += (pi0_boot[j] - min_pi0) * (pi0_boot[j] - min_pi0)
        min_mse = min(mse)
        argmin_mse = [i for i, mse_i in enumerate(mse) if mse_i == min_mse]
        pi0 = min(pi0[i] for i in argmin_mse)
        pi0 = min(pi0, 1)
    else:
        try:
            l = float(l)
        except ValueError:
            raise TypeError("lambda must be a number")
        if l < 0 or l >= 1:
            raise ValueError("lambda must be within [0,1)")
        pi0 = get_pi0(pv, [l])
        pi0 = min(pi0[0], 1)

    qvalues = m * [0]
    mincoeff = pi0 * pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in xrange(m-2, -1, -1):
        coeff = pi0*m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    return qvalues

def sparse_sub_dot(sparse_mtx1, sparse_mtx2, sub_size_x, sub_x, sub_size_y, sub_y):
  """ 
    multiply two submatrix to generate the dot product of 
    XY[sub_x*sub_size_x:(sub_x+1)*sub_size_x, sub_y*sub_size_y:(sub_y+1)*sub_size_y]
    sub_x, sub_x_start, sub_x_end: 0-based order of submatrix, its [start, end) indices
    sub_y, sub_y_start, sub_y_end: 0-based order of submatrix, its [start, end) indices
  """

  sub_x_start = sub_x * sub_size_x
  if(sparse_mtx1.shape[0]>(sub_x+1)*sub_size_x):
    sub_x_end=(sub_x+1)*sub_size_x
  else:
    sub_x_end=sparse_mtx1.shape[0]
  sub_y_start = sub_y * sub_size_y
  if(sparse_mtx2.shape[1]>(sub_y+1)*sub_size_y):
    sub_y_end=(sub_y+1)*sub_size_y
  else:
    sub_y_end=sparse_mtx2.shape[1]
  sub_X = sparse_mtx1.tocsr()[sub_x_start:sub_x_end,:]
  sub_Y = sparse_mtx2.tocsc()[:,sub_y_start:sub_y_end]
  sub_M = sub_X.dot(sub_Y)

  return sub_M

def sparse_std(sparse_mtx):
  a = sparse_mtx.copy(); a.data **= 2; return numpy.sqrt(a.mean()-sparse_mtx.mean()**2)

def is_symmetric(m):
  """Check if a sparse matrix is symmetric

  Parameters
  ----------
  m : array or sparse matrix
      A square matrix.

  Returns
  -------
  check : bool
      The check result.

  """
  if m.shape[0] != m.shape[1]:
    raise ValueError('m must be a square matrix')

  if not isinstance(m, coo_matrix):
    m = coo_matrix(m)

  r, c, v = m.row, m.col, m.data
  tril_no_diag = r > c
  triu_no_diag = c > r

  if triu_no_diag.sum() != tril_no_diag.sum():
    return False

  rl = r[tril_no_diag]
  cl = c[tril_no_diag]
  vl = v[tril_no_diag]
  ru = r[triu_no_diag]
  cu = c[triu_no_diag]
  vu = v[triu_no_diag]

  sortl = numpy.lexsort((cl, rl))
  sortu = numpy.lexsort((ru, cu))
  vl = vl[sortl]
  vu = vu[sortu]

  check = numpy.allclose(vl, vu)

  return check

def mad(data, axis=None):
  return numpy.median(numpy.absolute(data - numpy.median(data, axis)), axis)

def csr_median_mad(csr):
  csr_zeros = csr.shape[0]*csr.shape[1] - csr.nnz
  if csr_zeros >= 0.5*csr.shape[0]*csr.shape[1]:
    csr_median = 0
  else:
    csr_median_index = int(0.5*csr.shape[0]*csr.shape[1]) - csr_zeros
    csr_median = numpy.partition(csr.data, csr_median_index)[csr_median_index]
  deviate_zero = numpy.absolute(csr_median)
  deviate_data = numpy.absolute(csr.data - csr_median)
  deviate_data = numpy.insert(deviate_data, 0, deviate_zero)
  rank_zero = numpy.argwhere( numpy.argsort( deviate_data ) == 0 )[0,0]  #find rank of 0
  # d0, d1, d2, ...., dn
  # d1, d2, ..., d0(m), dm+1, ..., dn
  if rank_zero > int(0.5*csr.shape[0]*csr.shape[1]):
    csr_mad = numpy.sort(deviate_data[1:])[ int(0.5*csr.shape[0]*csr.shape[1]) ]
  elif rank_zero + csr_zeros < int(0.5*csr.shape[0]*csr.shape[1]):
    csr_mad = numpy.sort(deviate_data[1:])[ int(0.5*csr.shape[0]*csr.shape[1]) - csr_zeros + 1 ]
  else:
    csr_mad = deviate_zero
  return (csr_median, csr_mad)
  #scipy.sparse.random(m, n, density=0.01, format='coo', dtype=None, random_state=None, data_rvs=None)
  #csr_rand6=scipy.sparse.rand(60, 60, density=0.6, format='csr')

# A.indptr is an array, one for each row (+1 for the nnz):
def csr_row_set_nz_to_val(csr, row, value=0):
  """Set all nonzero elements (elements currently in the sparsity pattern)
  to the given value. Useful to set to 0 mostly.
  """
  if not isinstance(csr, scipy.sparse.csr_matrix):
    raise ValueError('Matrix given must be of CSR format.')
  csr.data[csr.indptr[row]:csr.indptr[row+1]] = value

# A.indptr is an array, one for each row (+1 for the nnz):
def csc_col_set_nz_to_val(csc, col, value=0):
  """Set all nonzero elements (elements currently in the sparsity pattern)
  to the given value. Useful to set to 0 mostly.
  """
  if not isinstance(csc, scipy.sparse.csc_matrix):
    raise ValueError('Matrix given must be of CSC format.')
  csc.data[csc.indptr[col]:csc.indptr[col+1]] = value

# Now you can just do:
#for row in indices:
#    csr_row_set_nz_to_val(A, row, 0)

# And to remove zeros from the sparsity pattern:
#A.eliminate_zeros()

def save_sparse_csr(filename,array):
  numpy.savez(filename, data = array.data, indices = array.indices,
                        indptr = array.indptr, shape = array.shape)

def load_sparse_csr(filename):
  loader = numpy.load(filename)
  return csr_matrix((  loader['data'], loader['indices'], loader['indptr']),
                       shape = loader['shape'])

def bed2df(bedfile, labels, chunksize=1024*1024):
  reader = pybedtools.BedTool(bedfile).to_dataframe(comment='#',names=labels,low_memory=False,sep="\t",chunksize=chunksize)
  i = 0
  dfs = []
  for chunk in reader:
    dfs.append(chunk)
    i += 1
    #NOTE:
    #               to avoid "pandas: Error tokenizing data. C error: out of memory"
  clust = pandas.concat( dfs )
  clust = clust.reindex()
  return clust

"""
  lib10x.getDistr(func, *args, **kwargs)
  ----------------------------------------------------
"""

def getDistr(key, distr):
  if key <= 0:
    return 1   #invalid key, not affecting ratio
  try:
    return distr[key]
  except:
    return distr[distr.last_valid_index()]

"""
  lib10x.cp2g(func, *args, **kwargs)
  ------------------------------------------------------------------------

  Convert from chromosome position to global index

    kwargs:
      chr=string: chromosome name
      pos=integer: position at chromosome
      A=dict: indexed accumulated chromosome length, exclusive

    output:
      integer 64: one number global index to encode both chromosome and positions

  We will need to optimize the speed of barrp data structure
"""

def cp2g(chr,pos,C):

  return C[chr]+pos-1

"""
  lib10x.g2cp(func, *args, **kwargs)
  ------------------------------------------------------------------------

  Convert global index back to chromosome and position

    kwargs:
      gin=integer 64: a 64-bit integer to encode global genomic coordinates
      A=dict: indexed accumulated chromosome length, exclusive
      acclen=list: accumulated chromosome length
      chrname=list: ordered list of chromosome names

    output:
      (string,int): chrom, position

  We will need to optimize the speed of barrp data structure
"""

def g2cp(gin,C,acclen_C,chrname):

  chr=chrname[bisect.bisect(acclen_C,gin)-1]
  """ if gin at acclen, give right, logN """

  return (chr,gin-C[chr]+1)

"""
  lib10x.int2bool(func, num)
  ---------------------------------------------------------------

  Convert a 0 to 7 integer to a boolean list, e.g.:
    1 -> [T, F, F]
    7 -> [T, T, T]

"""

def int2bool(num):
  return [bool(num & (1<<n)) for n in range(3)]

"""
  lib10x.load_refdata(func, *args, **kwargs)
  ------------------------------------------------------------------------
"""

def load_refdata(bam):

  samfile = pysam.AlignmentFile(bam, "rb")
  refdata=samfile.header['SQ']
  # A[chr] contains accumulated size up to 'chr'; C[chr] contains the size of the chromosome
  # cp2g: g = A[chr]+pos
  # g2cp: chr = max_c : A[c]<g; pos=g-A[chr]
  chrname=map(lambda t : t['SN'], refdata)
  chrlen=map(lambda t : t['LN'], refdata)
  acclen=list(numpy.cumsum(chrlen))
  A=collections.OrderedDict(zip(chrname,acclen))
  B=collections.OrderedDict(zip(chrname,acclen[1:] + [2*acclen[-1]]))
  refdata={'chrname':chrname,'chrlen':chrlen,'acclen':acclen,'A':A,'B':B}

  return refdata

"""
  lib10x.load_quad(func, *args, **kwargs)
  ------------------------------------------------------------------------

  Load barcoded reads from an input quad.

    kwargs:
      bam=bamfile: 10x bamfile
      bed=BedTool: interval region
      bar=barcode: barcode set

    output:
      barrp: {bar:{'c1':[br1,br2,...,brn],'p1':[br1,br2,...,brn],'c2':[br1,br2,...,brn],'p2':[br1,br2,...,brn]},...}

  We will need to optimize the speed of barrp data structure
"""

def load_quad(bam,bed,bar,i):

  print >>sys.stderr, "loading ", i, "-th region"
  chrom,start,end=bed
  samfile = pysam.AlignmentFile(bam, "rb")
  N_all = 0; N_qc=0; N_bar=0;
  barrp = dict(zip(bar,[{'c1':[],'c2':[],'p1':[],'p2':[]} for _ in range(len(bar))]))
  #NOTE: *len(bar) is not working here because * only replicates the reference

  for rd in samfile.fetch(chrom,start,end):
    #if N_all % 100000 == 0 and N_all != 0: print >>sys.stderr, N_all, "reads processed"
    #if N_all >= maxread: break
    N_all += 1

    if rd.is_unmapped or rd.is_qcfail or rd.mapping_quality<30 or rd.mate_is_unmapped or rd.is_secondary: #N_all
      #throw unmapped and low quality read and qc failed read read pairs
      continue
    bar = dict(rd.tags).get('BX')
    if not bar:
      #throw Barcode unresolved read pairs
      continue

    N_qc += 1
    bar = bar.split('-')[-1] #BX is the barcode field

    #we then gather read pair data only from read1 trusting mrnm and mpos
    if rd.is_read1:
      if not barrp.has_key(bar): #barcode is new
        continue #ignore anything not in requested barcode
      else: #barcode exists
        N_bar += 1
        barrp[bar]['c1'].append(rd.rname)
        barrp[bar]['p1'].append(rd.pos)
        barrp[bar]['c2'].append(rd.mrnm)
        barrp[bar]['p2'].append(rd.mpos)

  return barrp

"""
  lib10x.compute_uv(func, *args, **kwargs)
  ------------------------------------------------------------------------

  Compute likelihood from four-quad input barrp (u-,u+), (v-,v+)

    Assumptions:
      region u preceeds region v in genome
      loaded barcode sets in u and v are exactly the same

    kwargs:
      u=reg: region u
      v=reg: region v
      urp=barrp: a pair of barrp data from load_quad ordered (u-,u+)
      vrp=barrp: a pair of barrp data from load_quad ordered (v-,v+)
      alt=alternative: C i.e. (u-,v-) and (u+,v+)
                       X i.e. (u-,v+) and (v-,u+)
      mode=likelihood: integer 1-7 for at least one of (H,G,T)
                       1st bit H i.e. likelihood from fragment length density
                       2nd bit G i.e. likelihood from intra-fragment spacing density
                       3rd bit T i.e. likelihood from barcode frequency density

    output:
      a log likelihood ratio

  We will need to optimize the speed of mathematical calculation
"""


def compute_uv(u,v,urp,vrp,DistrH,DistrG,DistrT,alt='C',mode=7,pi=0):

  #assume w-x, y-z and w<y, w<x, y<z (by coordinates)
  #w always the smallest of four quads, i.e. u-
  #x is either u+ as for N or v- as for C or v+ as for X
  #y is either v- for N or u+ for C and X
  #z is either v+ as for N and C or v- as for X
  #u is middle point of u- and u+
  #v is middle point of v- and v+

  print >>sys.stderr, "compute_uv", "u=", u, "v=", v
  w1 = urp[0]; x1 = urp[1]; y1=vrp[0]; z1=vrp[1]
  w1c = -u; x1c = u; y1c = -v; z1c = v
  if alt == 'C':
    w2 = urp[0]; x2 = vrp[0]; y2=urp[1]; z2=vrp[1]
    w2c = -u; x2c = -v; y2c = u; z2c = v
  elif alt == 'X':
    w2 = urp[0]; x2 = vrp[1]; y2=urp[1]; z2=vrp[0]
    w2c = -u; x2c = v; y2c = u; z2c = -v
  else:
    quit()
    #for b in urp[0].keys(): print b,":",urp[0][b]['p1']
    #for b in urp[1].keys(): print b,":",urp[1][b]['p1']
    #for b in vrp[0].keys(): print b,":",vrp[0][b]['p1']
    #for b in vrp[1].keys(): print b,":",vrp[1][b]['p1']

  loglikelihoodratio = 0.0

  use_lr = dict(zip(['H','G','T'],int2bool(mode))) # {H:Y/N, G:Y/N, T:Y/N}

  nw1x1=0; ny1z1=0; nw2x2=0; ny2z2=0 #calculate total number of molecules given connection
  for b in w1.keys():
    #print len(w1[b]['p1']), len(x1[b]['p1']), len(y1[b]['p1']), len(z1[b]['p1'])
    if len(w1[b]['p1'])>0 or len(x1[b]['p1'])>0:  nw1x1 += 1;
    if len(y1[b]['p1'])>0 or len(z1[b]['p1'])>0:  ny1z1 += 1;
    if len(w2[b]['p1'])>0 or len(x2[b]['p1'])>0:  nw2x2 += 1;
    if len(y2[b]['p1'])>0 or len(z2[b]['p1'])>0:  ny2z2 += 1;
  #print >>sys.stderr, "total fragments", "nw1x1=", nw1x1, "ny1z1=", ny1z1, "nw2x2=", nw2x2, "ny2z2=", ny2z2
  #quit()

  for b in w1.keys(): #iterate through barcodes
    print "b=", b
    print "r=", w1[b]['p1'],x1[b]['p1'],y1[b]['p1'],z1[b]['p1'],w2[b]['p1'],x2[b]['p1'],y2[b]['p1'],z2[b]['p1']
    if use_lr['H']:
      lrh = float(compute_lrh(w1[b]['p1'],x1[b]['p1'],y1[b]['p1'],z1[b]['p1'],w2[b]['p1'],x2[b]['p1'],y2[b]['p1'],z2[b]['p1'],DistrH))
      print "lrh=", lrh
    else:
      lrh = 1
    if use_lr['G']:
      lrg = float(compute_lrg(w1[b]['p1'],x1[b]['p1'],y1[b]['p1'],z1[b]['p1'],w2[b]['p1'],x2[b]['p1'],y2[b]['p1'],z2[b]['p1'],w1c,x1c,y1c,z1c,w2c,x2c,y2c,z2c,DistrG))
      print "lrg=", lrg
    else:
      lrg = 1
    if use_lr['T']:
      lrt = float(compute_lrt(w1[b]['p1'],x1[b]['p1'],y1[b]['p1'],z1[b]['p1'],w2[b]['p1'],x2[b]['p1'],y2[b]['p1'],z2[b]['p1'],nw1x1,ny1z1,nw2x2,ny2z2,b,DistrT))
      print "lrt=", lrt
    else:
      lrt = 1
    print "lrh*lrg=", lrh*lrg
    loglikelihoodratio += numpy.log((1-pi)+pi*lrh*lrg*lrt)

  return loglikelihoodratio


"""
  lib10x.compute_lrh(func, *args, **kwargs)
  ------------------------------------------------------------------------

  Compute fragment length density based likelihood for change from null:(w1-x1,y1-z1) to alt:(w2-x2,y2-z2)

    Assumptions:

    kwargs:
      w1b,x1b,y1b,z1b: read1 positions of barcoded read sets under null
      w2b,x2b,y2b,z2b: read1 positions of barcoded read sets under alternative
      DistH: fragment length density

    output:
      a likelihood ratio

  We will need to optimize the speed of mathematical calculation
"""

def compute_lrh(w1b,x1b,y1b,z1b,w2b,x2b,y2b,z2b,DistrH):
  return getDistr(len(w2b)+len(x2b),DistrH)*getDistr(len(y2b)+len(z2b),DistrH)/(getDistr(len(w1b)+len(x1b),DistrH)*getDistr(len(y1b)+len(z1b),DistrH))


"""
  lib10x.compute_lrg(func, *args, **kwargs)
  ------------------------------------------------------------------------

  Compute read spacing density based likelihood for change from null:(w1-x1,y1-z1) to alt:(w2-x2,y2-z2)

    Assumptions:

    kwargs:
      w1b,x1b,y1b,z1b: read1 positions of barcoded read sets under null
      w2b,x2b,y2b,z2b: read1 positions of barcoded read sets under alternative
      w1c,x1c,y1c,z1c: cut point next to the null read set, sign '-' if cut is to the right, '+' if cut to the left
      w2c,x2c,y2c,z2c: cut point next to the alt read set, sign '-' if cut is to the right, '+' if cut to the left
      DistG: read spacing density

    output:
      a likelihood ratio

  We will need to optimize the speed of mathematical calculation
"""

def compute_lrg(w1b,x1b,y1b,z1b,w2b,x2b,y2b,z2b,w1c,x1c,y1c,z1c,w2c,x2c,y2c,z2c,DistrG):
  if(len(w1b)>1 and len(x1b)>1):
    dw1x1 = int(0.5 * (numpy.sign(w1c)+1)*(min(w1b)-abs(w1c)) - 0.5 * (numpy.sign(w1c)-1)*(abs(w1c)-max(w1b)) \
            + 0.5 * (numpy.sign(x1c)+1)*(min(x1b)-abs(x1c)) - 0.5 * (numpy.sign(x1c)-1)*(abs(x1c)-max(x1b)))
  else:
    dw1x1 = 0

  if(len(y1b)>1 and len(z1b)>1):
    dy1z1 = int(0.5 * (numpy.sign(y1c)+1)*(min(y1b)-abs(y1c)) - 0.5 * (numpy.sign(y1c)-1)*(abs(y1c)-max(y1b)) \
            + 0.5 * (numpy.sign(z1c)+1)*(min(z1b)-abs(z1c)) - 0.5 * (numpy.sign(z1c)-1)*(abs(z1c)-max(z1b)))
  else:
    dy1z1 = 0

  if(len(w2b)>1 and len(x2b)>1):
    dw2x2 = int(0.5 * (numpy.sign(w2c)+1)*(min(w2b)-abs(w2c)) - 0.5 * (numpy.sign(w2c)-1)*(abs(w2c)-max(w2b)) \
            + 0.5 * (numpy.sign(x2c)+1)*(min(x2b)-abs(x2c)) - 0.5 * (numpy.sign(x2c)-1)*(abs(x2c)-max(x2b)))
  else:
    dw2x2 = 0

  if(len(y2b)>1 and len(z2b)>1):
    dy2z2 = int(0.5 * (numpy.sign(y2c)+1)*(min(y2b)-abs(y2c)) - 0.5 * (numpy.sign(y2c)-1)*(abs(y2c)-max(y2b)) \
            + 0.5 * (numpy.sign(z2c)+1)*(min(z2b)-abs(z2c)) - 0.5 * (numpy.sign(z2c)-1)*(abs(z2c)-max(z2b)))
  else:
    dy2z2 = 0

  return getDistr(dw2x2,DistrG)*getDistr(dy2z2,DistrG)/(getDistr(dw1x1,DistrG)*getDistr(dy1z1,DistrG))

"""
  lib10x.compute_lrt(func, *args, **kwargs)
  ------------------------------------------------------------------------

  Compute barcode fragment density based likelihood for change from null:(w1-x1,y1-z1) to alt:(w2-x2,y2-z2)

    Assumptions:

    kwargs:
      nw1x1, ny1z1: number of fragments under null
      nw2x2, ny2z2: number of fragments under alt
      b: a specific barcode
      DistT: barcode fragment density

    output:
      a likelihood ratio

  We will need to optimize the speed of mathematical calculation
"""

def compute_lrt(w1b,x1b,y1b,z1b,w2b,x2b,y2b,z2b,nw1x1,ny1z1,nw2x2,ny2z2,b,DistrT):
  pseudo=0.00000000001
  lw2x2 = nw2x2*DistrT[b]+pseudo; ly2z2 = ny2z2*DistrT[b]+pseudo; lw1x1 = nw1x1*DistrT[b]+pseudo; ly1z1 = ny1z1*DistrT[b]+pseudo;
  return poisson.pmf((len(w2b)+len(x2b))>0,lw2x2)*poisson.pmf((len(y2b)+len(z2b))>0,ly2z2)/poisson.pmf((len(w1b)+len(x1b))>0,lw1x1)/poisson.pmf((len(y1b)+len(z1b))>0,ly1z1)
