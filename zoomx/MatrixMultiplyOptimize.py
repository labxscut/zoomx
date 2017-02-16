#!/usr/bin/env python

import numpy
import scipy.sparse
import time
import sys

#60K regs, 60K bars, 33M frags
X = scipy.sparse.rand(60000, 60000, density=0.02, format='csr')
Y = X.transpose(copy=True)
print >>sys.stderr, "multiplying matrix X and Y..."
gtk1_time = time.time()
X
Y
M = X.dot(Y)
gtk2_time = time.time()
print >>sys.stderr, gtk2_time-gtk1_time, "s done"
