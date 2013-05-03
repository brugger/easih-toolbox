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

fasta_ref = pysam.Fastafile("/data/refs/HIV/K03455.fasta");


MIN_PERCENTAGE = 1

def codon2AA(codon):

    codon2AA = { 'TTT': 'F', 'TTC': 'F',
                 'TTA': 'L', 'TTG': 'L',
                 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                 'TAT': 'Y', 'TAC': 'Y', 
                 'TAA': '*', 'TAG': '*', 'TGA': '*',
                 'TGT': 'C', 'TGC': 'C',
                 'TGG': 'W',
                 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                 'CAT': 'H', 'CAC': 'H',
                 'CAA': 'Q', 'CAG': 'Q',
                 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
                 'ATG': 'M',
                 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                 'AAT': 'N', 'AAC': 'N',
                 'AAA': 'K', 'AAG': 'K',
                 'AGT': 'S', 'AGC': 'S',
                 'AGA': 'R', 'AGG': 'R',
                 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                 'GAT': 'D', 'GAC': 'D',
                 'GAA': 'E', 'GAG': 'E',
                 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G' }
                 
    if (len(codon) != 3):
        return "x"


    if ( codon not in codon2AA ):
#        print "Cannot tranlate " + codon + " codon to an AA";
        return "X"

    return codon2AA[ codon ] 



def codon_count2AA_count( codon_counts ):
    AA_counts = dict();

#    pp.pprint( codon_counts )
    

    total_AAs = 0

    for codon in (codon_counts.keys()):
        if (len(codon) != 3):
            continue

        AA = codon2AA( codon )

        if ( AA == "X"):
            continue

        if ( AA not in AA_counts) :
            AA_counts[ AA ] = 0;

        AA_counts[ AA ] += codon_counts[ codon ];
        total_AAs += codon_counts[ codon ]


#    pp.pprint( AA_counts )

    AA_percents = dict()

    for AA in (AA_counts.keys()):
        AA_percent = float(100)*float(AA_counts[ AA ])/float( total_AAs )
        if ( AA_percent >= MIN_PERCENTAGE ):
            AA_percents[ AA ] = "%.2f" % ( AA_percent )

#    pp.pprint( AA_percents )

    return AA_percents

#############################
def update_counts(codon_counts, codon):

    if (codon not in codon_counts):
        codon_counts[ codon ] = 0

    codon_counts[ codon ] += 1

    return codon_counts;


Stamford_format = dict();


gene_list = [['Pol', 2253, 2546],['RT', 2550, 4227]]
#gene_list = [['Test', 2550, 4227]]

for gene in gene_list:
    (gene_name, ORF_start, ORF_end) = (gene[0], gene[1], gene[2]);

    Stamford_format[gene_name] = [];

#    print "\t".join([gene_name, str(ORF_start), str(ORF_end)])
    print gene_name

    iter = samfile.pileup("K03455", ORF_start -1 , ORF_end - 1, max_depth=80000)

    skip_bases = 0
    exit_counter = 50

    READING_FRAME_OFFSET = 0

    INDEL_SKIPPING = 0

    for x in iter:
 

        if ( x.pos < ORF_start - 1):
            continue

        if ( x.pos > ORF_end ):
            break


        if ( INDEL_SKIPPING ):
                #        print "skipping an indel base: " + str( INDEL_SKIPPING )
            INDEL_SKIPPING += 1
            continue


        indel_counter = dict()

        for read in x.pileups:
            
            if ( read.indel not in indel_counter) :
                indel_counter[ read.indel ] = 0;
                
                indel_counter[ read.indel ] += 1;
                
                
        major_mutation = sorted(indel_counter, key=indel_counter.get, reverse=True);

        if ( major_mutation[0] <0 ):
            INDEL_SKIPPING = major_mutation[0]

        READING_FRAME_OFFSET += major_mutation[0]

        if ( skip_bases ):
            skip_bases -= 1
            continue
   
        skip_bases = 2
    
        passed_bases = 0

        codon_counts = dict()

        for read in x.pileups:
            if (read.alignment.is_unmapped or read.alignment.is_duplicate or read.is_del):
                continue

 
            if ( read.indel > 0):
                insertion = read.alignment.seq[ read.qpos:read.qpos+ read.indel + 1]
                codon_counts = update_counts(codon_counts, insertion)
             
             

            elif ( read.indel < 0):
                deletion = ''
                deletion = read.alignment.seq[ read.qpos] + "-" * abs(read.indel)

                codon_counts = update_counts(codon_counts, deletion)
            
            
            else:
                alt_base = read.alignment.seq[ read.qpos:read.qpos+3];
                if (len( alt_base) != 3):
                    continue
                codon_counts = update_counts(codon_counts, alt_base)


#    pp.pprint( codon_counts )
        

        good_AA_counts = codon_count2AA_count( codon_counts )

        reference_codon = fasta_ref.fetch(str(samfile.getrname(x.tid)), x.pos + READING_FRAME_OFFSET, x.pos+3 + READING_FRAME_OFFSET)
        reference_AA    = codon2AA( reference_codon)

        codon_number    = (1+ (x.pos + 1 - ORF_start)/3)

        if ( len(good_AA_counts.keys()) >= 2):

            frequencies = []
            for AA in (sorted(good_AA_counts.keys())):
                if ( reference_AA == AA):
                    continue;

                frequencies.append( AA+":"+good_AA_counts[ AA ]+"%")

                Stamford_format[gene_name].append("%s%d%s" % (reference_AA, codon_number, AA))
        
#        print "reading frame offset: " + str(READING_FRAME_OFFSET)

            print str(x.pos+1) + "\t" + reference_AA +"%d\t" % codon_number + "\t".join( frequencies )


        else:

            if (reference_AA in good_AA_counts):
                continue

            frequencies = []
            for AA in (sorted(good_AA_counts.keys())):
                frequencies.append( AA+":"+good_AA_counts[ AA ]+"%")
        
#        print "reading frame offset: " + str(READING_FRAME_OFFSET)

                Stamford_format[gene_name].append("%s%d%s" % (reference_AA, codon_number, AA))

            print str(x.pos+1) + "\t" + reference_AA +"%d\t" % codon_number + "\t".join( frequencies )


        if (0):
            exit_counter -= 1
            if ( not exit_counter):
                exit()
            

    print 

print 
print    

for gene in (sorted(Stamford_format.keys())):
    print gene
    print ", ".join(Stamford_format[ gene ]);

