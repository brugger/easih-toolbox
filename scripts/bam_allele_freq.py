#!/usr/bin/python
import sys
import getopt
import pprint
pp = pprint.PrettyPrinter(indent=4)

sys.path.append('/home/kb468/src/pysam-0.6/build/lib.linux-x86_64-2.6/')
sys.path.append('/scratch/kb468/pysam-0.6/build/lib.linux-x86_64-2.6/')

import re

import pysam


bamfile = sys.argv[1]
samfile = pysam.Samfile( bamfile, "rb" )
#samfile = pysam.Samfile( "CP0007/C070003-1.bam", "rb" )


chromosome = str(sys.argv[2])
position   = int(sys.argv[3])
position -= 1

#############################
def update_counts(base_counts, base):

    if (base not in base_counts):
        base_counts[ base ] = 0

    base_counts[ base ] += 1

    return base_counts;



iter = samfile.pileup(max_depth=20000)

iter = samfile.pileup(chromosome, position,position+1,max_depth=20000);

MIN_BASEQ        =  0


for x in iter:
 
    if (position != x.pos ):
        continue
    
    base_counts = dict()

#    print 'coverage at base %s = %s' %(x.pos , x.n)

    for read in x.pileups:
        if (read.alignment.is_unmapped or read.alignment.is_duplicate or read.is_del):
            continue

        if ( read.indel > 0):
            insertion = read.alignment.seq[ read.qpos:read.qpos+ read.indel + 1]
            base_counts = update_counts(base_counts, insertion)
             
        elif ( read.indel < 0):
            deletion = ''
            deletion = read.alignment.seq[ read.qpos] + "-" * abs(read.indel)
            deletion = deletion.lower()
            deletion = deletion.title()

            base_counts = update_counts(base_counts, deletion)
            
        else:
            base_qual = ord( read.alignment.qual[ read.qpos] ) - 33
            if ( base_qual < MIN_BASEQ):
                continue
            alt_base = read.alignment.seq[ read.qpos ];
            base_counts = update_counts(base_counts, alt_base)


#    pp.pprint( base_counts )

    results = []
    results.append( "%s:%d" % (chromosome,position+1))

    for allele in base_counts.keys():
        results.append('%s:%d' %( allele, base_counts[ allele]))


    print "\t".join( results )

samfile.close()

