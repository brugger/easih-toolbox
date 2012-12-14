#!/usr/bin/python
import sys
import getopt
import pprint
pp = pprint.PrettyPrinter(indent=4)

import scipy
from scipy import stats

import pysam

bamfile = sys.argv[1]
samfile = pysam.Samfile( bamfile, "rb" )
#samfile = pysam.Samfile( "CP0007/C070003-1.bam", "rb" )

chromosome = str(sys.argv[2])
position   = int(sys.argv[3])
position   -= 1

MIN_MAPQ         =  0
MIN_BASEQ        =  0


def update_counts(base_counts, alt, qual, strand, pos):

    if ( alt not in base_counts):
        base_counts[ alt ] = dict()
        base_counts[ alt ][ 'count' ] = 0
        base_counts[ alt ][ 'quals' ] = 0
        base_counts[ alt ][ 0 ] = dict()
        base_counts[ alt ][ 1 ] = dict()

    base_counts[ alt ][ 'count' ] += 1
    base_counts[ alt ][ 'quals' ] += qual

    if (pos not in base_counts[ alt ][ strand ]):
        base_counts[ alt ][ strand ][ pos ] = 1
    else:
        base_counts[ alt ][ strand ][ pos ] = +1

    return base_counts;



fasta_ref = pysam.Fastafile("/data/refs/human_1kg/human_g1k_v37.fasta");
iter = samfile.pileup(chromosome, position,position+1,max_depth=20000);

for x in iter: 
#    print str(x)

    if ( position != x.pos):
        continue

    ref_base = fasta_ref.fetch(str(samfile.getrname(x.tid)), x.pos, x.pos+1 )


    passed_bases = 0

    base_counts = dict()

    for read in x.pileups:

        if (read.alignment.is_unmapped or read.alignment.is_duplicate or read.is_del):
            continue

        if (read.alignment.mapq <= MIN_MAPQ ):
            continue

        if ( read.indel > 0):
            insertion = read.alignment.seq[ read.qpos:read.qpos+ read.indel + 1]
            base_qualities = 0

            for base_qual in read.alignment.qual[ read.qpos:read.qpos+ read.indel + 1]:
                base_qualities += ord( base_qual ) - 33

            base_qual = base_qualities/(read.indel + 1)

            if ( base_qual < MIN_BASEQ):
                continue

            base_counts = update_counts(base_counts, insertion, 
                                        base_qual, read.alignment.is_reverse, read.alignment.aend)

#            pp.pprint( base_counts )
        elif ( read.indel < 0):
            deletion = fasta_ref.fetch(str(samfile.getrname(x.tid)), x.pos, x.pos+abs(read.indel)+1 )

            base_qual = ord( read.alignment.qual[ read.qpos] ) - 33

            if ( base_qual < MIN_BASEQ):
                continue

            base_counts = update_counts(base_counts, deletion, 
                                        base_qual, read.alignment.is_reverse, read.alignment.aend)
        else:
            base_qual = ord( read.alignment.qual[ read.qpos] ) - 33

            if ( base_qual < MIN_BASEQ):
                continue

            alt_base = read.alignment.seq[ read.qpos];
            base_counts = update_counts(base_counts, alt_base, 
                                        base_qual, read.alignment.is_reverse, read.alignment.aend)

#    pp.pprint( base_counts )


    print "Position: " + samfile.getrname(x.tid) +":"+ str(x.pos + 1)

    for allele in base_counts.keys():
        print "\t".join([ allele,
                          "reads: " + str(base_counts[ allele ][ 'count' ]),
                          "avg qual: %.2f" % (float(base_counts[ allele ][ 'quals' ])/ float(base_counts[ allele ][ 'count' ])),
                          "forward starts: " +  str(len(base_counts[ allele ][ 0 ].keys())),
                          "reverse starts: " +  str(len(base_counts[ allele ][ 1 ].keys())),
                          
                          "."
                          ])

        


    first_allele = 'N'
    first_count  = -1

    second_allele = 'N'
    second_count  = -1

#    if ( x.pos == 135779171 ):
#        pp.pprint( base_counts )

    for allele in base_counts.keys():
        if (base_counts[ allele ][ 'count'] > first_count ):
            
            second_allele = first_allele 
            second_count  = first_count
            
            first_allele = allele
            first_count   = base_counts[ allele ][ 'count']
            
        elif (base_counts[ allele ][ 'count'] > second_count):
            second_allele = allele
            second_count = base_counts[ allele ][ 'count']
            
            
    first_second_ratio = float( second_count)/float(first_count + second_count)*100


    if ( first_allele == 'N' or second_allele == 'N'):
        continue

    zvalue,pvalue2 = stats.ranksums(base_counts[ first_allele ][ 0 ].keys()+ base_counts[ first_allele ][ 1 ].keys(),
                                    base_counts[ second_allele ][ 0 ].keys()+ base_counts[ second_allele ][ 1 ].keys())



            # Am using fishers test to see if there is a strand bias for the SNP
    oddsratio, pvalue = stats.fisher_exact([[len(base_counts[ first_allele  ][ 0 ].keys()), 
                                             len(base_counts[ first_allele  ][ 1 ].keys())], 
                                            [len(base_counts[ second_allele ][ 0 ].keys()), 
                                             len(base_counts[ second_allele ][ 1 ].keys())]])

    
    print "Major allele: " + first_allele
    print "Minor allele: " + second_allele
    print "Allele ratio: %.2f %%" % first_second_ratio
    
    print "strand bias (pvalue): %.2f (> 0.05 no bias)" % pvalue
    print "read start bias (pvalue): %.2f (> 0.05 no bias)" % pvalue2



            

