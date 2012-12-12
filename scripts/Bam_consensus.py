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
#samfile = pysam.Samfile( "../miseq/A280050.bam", "rb" )

#samfile = pysam.Samfile( "/data/CP/CP0007/C060004.bam", "rb" )


MIN_ALLELE_PERC = 25



#############################
def First_Char(base):
    """Return only first character of a multi-char base"""
    
    if(base):
        
            bs = re.search("^(.)", base)
            if(bs):
                b = bs.group(1)
    return(b)



#############################
def Find_IUPAC(fb,sb,tb = ''):
    """Return an IUPAC code given a first and second base"""


#    print "\t".join([fb,sb,tb])


    max_length =  max(len(fb), len(sb), len(tb))

    IUPAC_Codes = {'A': 'A',
                   'C': 'C',
                   'G': 'G',
                   'T': 'T',
                   'AT': 'W',
                   'CG': 'S',
                   'AC': 'M',
                   'GT': 'K',
                   'AG': 'R',
                   'CT': 'Y',
                   'CGT': 'B',
                   'AGT': 'D', 
                   'ACT': 'H',
                   'ACG': 'V'
                   }

    IUPAC_consensus = ""

    for i in range(0, max_length):
        bases = [];

#        print str(len(fb)) , " >= " ,str(i)

        if (len(fb) > i):
            bases.append( fb[i] )

        if (len(sb) > i):
            bases.append( sb[i] )

        if (len(tb) > i):
            bases.append( tb[i] )

        fst = ''.join( sorted(set(bases)))

        IUPAC_consensus += IUPAC_Codes[ fst ]


#    print IUPAC_consensus

    return(IUPAC_consensus)


#############################
def update_counts(base_counts, base):

    if (base not in base_counts):
        base_counts[ base ] = dict()

        base_counts[ base ][ 'count' ] = 0

    base_counts[ base ][ 'count' ] += 1

    return base_counts;



iter = samfile.pileup(max_depth=80000)
#iter = samfile.pileup("NC_001802", 1621, 1621 + 10, max_depth=80000)

consensus = []

hyphenlen = 0 ###

for x in iter:
 
#############################
###Fix indel###
    
    if(hyphenlen):
        hyphenlen -= hyphenlen
        continue

#############################

    
    passed_bases = 0

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
            alt_base = read.alignment.seq[ read.qpos];
            base_counts = update_counts(base_counts, alt_base)


#    pp.pprint( base_counts )


    #continue

    
    base_perc = dict()
    first_base    = 'N'
    first_count  = -1

    second_base  = 'N'
    second_count  = -1

    third_base  = 'N'
    third_count  = -1

   
    for base in sorted (base_counts, key = base_counts.get, reverse = True):
 #       print "%s, %s" %(base , base_counts[ base ])
        if (base_counts[ base ][ 'count'] > first_count ):

            second_base   = first_base 
            second_count = first_count

            first_base = base
            first_count = base_counts[ base ][ 'count']
            
        elif (base_counts[ base ][ 'count'] > second_count):

            third_base = second_base
            third_count = second_count
            
            second_base = base
            second_count = base_counts[ base ][ 'count']
            
        elif (base_counts[ base ][ 'count'] > third_count):
            third_base = base
            third_count = base_counts[ base ][ 'count']


###############################            
### Indel fix ###            
    if(len(first_base) > 1):
        hyphen = re.search("^([ATGCNatgcn])([-]+)", first_base)
            
        if(hyphen):
            anchor = hyphen.group(1)
            first_base = anchor
            hyphenlen = len(hyphen.group(2))

###############################            


            

    ### If coverage is more than 4000
    if( x.n > 4000 ):

#        print "First count, second count, third count: %s, %s, %s" %(first_count, second_count, third_count)


        first_rest_perc = float(first_count) / (float(third_count) + float(second_count) + float(first_count))*100
#        first_rest_perc = round(first_rest_perc, 2)
        first_rest_perc = round(first_rest_perc)     
#        print "first_rest_perc: %s" %(first_rest_perc)

        second_rest_perc = float(second_count) / (float(third_count) + float(second_count) + float(first_count))*100
#        second_rest_perc = round(second_rest_perc, 2)
        second_rest_perc = round(second_rest_perc)
#        print "second_rest_perc: %s" %(second_rest_perc)

        third_rest_perc = float(third_count) / (float(third_count) + float(second_count) + float(first_count))*100
#        third_rest_perc = round(third_rest_perc, 2)
        third_rest_perc = round(third_rest_perc)
#        print "third_rest_perc: %s" %(third_rest_perc)
        
#        print "\t".join( [str(x.pos), str(first_rest_perc), str(second_rest_perc), str(third_rest_perc), first_base, second_base, third_base])


        ### Get consensus - based on base percentage
        if( third_rest_perc >= MIN_ALLELE_PERC ):

            new_base = Find_IUPAC(first_base,second_base,third_base)
            consensus.append(new_base)

                
        elif( second_rest_perc >= MIN_ALLELE_PERC ):

            new_base = Find_IUPAC(first_base,second_base)
            consensus.append(new_base)
        else:
            consensus.append(first_base)
            
            
#### get input filename without the path
b = re.search(".*\/?(.*)$", bamfile)
#print b
if ( b.group(1) ):
    bamfilename = b.group(1)
else:
    bamfilename = bamfile

print ">%s_Consensus" %(bamfilename)

#print ">IUPAC_Consensus"

dashlen = 0
i = 0

#### print consensus sequence
for j in range(len(consensus)):

    if j > i:
        i = j + dashlen
        dashlen = 0
        anchor = ''

          
        if(len(consensus[i]) > 1):
            dash = re.search("^([ATGCNatgcn])([-]+)", consensus[i])
                        
            if(dash):
                anchor = dash.group(1)
                dashlen = len(dash.group(2))
                sys.stdout.write(anchor)

            else:
                sys.stdout.write(consensus[i])

        else:
            sys.stdout.write(consensus[i])


print        


samfile.close()


###############################
    #IUPAC code 
    
    ## only 1 base represented at each position
    # ATGC

 #   if ( first_second_perc >= 75 ):
        
 #       consensus = first_base

    ## 2 bases represented at each position
    # W: A & T
    # S: G & C
    # M: A & C
    # K: G & T
    # R: A & G
    # Y: C & T


    ## 3 bases represented at each position
    # B: C, G & T
    # D: A, G & T
    # H: A, C & T
    # V: A, C & G

    ## 4 bases represented at each position
    # N or - : Any base (not a gap)
    

    
