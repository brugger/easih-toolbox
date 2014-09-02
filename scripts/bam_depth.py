#!/usr/bin/python
import sys
import getopt
import pprint
pp = pprint.PrettyPrinter(indent=4)

import scipy
from scipy import stats

sys.path.append('/software/congenica/packages/pysam-0.7.5/build/lib.linux-x86_64-2.7/')

import pysam

bamfile = sys.argv[1]
samfile = pysam.Samfile( bamfile, "rb" )

chromosome = str(sys.argv[2])
position   = int(sys.argv[3])
position   -= 1

MIN_MAPQ         =  0
MIN_BASEQ        =  0
iter = samfile.pileup(chromosome, position,position+1,max_depth=20000);

printed_depth = 0;

for x in iter: 

    if ( position != x.pos):
        continue

    depth = 0

    for read in x.pileups:

        if (read.alignment.is_unmapped or read.alignment.is_duplicate or read.is_del):
            continue

        if (read.alignment.mapq <= 0 ):
            continue

        depth += 1

    print "\t".join([samfile.getrname(x.tid), str(x.pos + 1), str(x.pos + 1), str(depth)])


