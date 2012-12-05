#!/usr/bin/python
import sys
import getopt

sys.path.append('/scratch/kb468/pysam-0.6/build/lib.linux-x86_64-2.6/')
import pysam

letters = 'b:B:' # the : means an argument needs to be passed after the letter
opts, extraparams = getopt.getopt(sys.argv[1:], letters) 

baits    = '/data/B20/NBG_V3_200.bed'
baits    = '/data/B20/NBG_V3.bed'
#baits    = '/scratch/kb468/pysam-0.6/test.bed'
#baits    = '/scratch/kb468/pysam-0.6/test2.bed'
bam_file = 'missing'

for o,p in opts:
  if o in ['-b']:
     bam_file = p
  elif o in ['-B']:
     baits = p

#print 'Baits: ', baits
#print 'Bam:   ', bam_file

verbose = 0

samfile = pysam.Samfile( bam_file, "rb" )

bedfile = open(baits, "r")
bedline = bedfile.readline().rstrip()
(chr, start, end) = bedline.split("\t")
if ( verbose ):
  print "|"+chr+"|\t|"+start+"|\t|"+end+"|"


mq0  = 0
mq20 = 0
mq50 = 0

empty_bedfile = 0;

for read in samfile.fetch():

  if ( read.is_unmapped ):
    continue

  read_chr   = samfile.getrname(read.tid)
  read_start = read.aend - read.alen  + 1
  read_end   = read_start + read.alen + 1

  if ( read_chr == "X"):
    read_chr = 23 

  if ( read_chr == "Y"):
    read_chr = 24 


  # So far only handling chromosome 1-22+X&Y(23 and 24)
  if ( not str(read_chr).isdigit()):
    if ( verbose >= 2):
      print "Not a real chromosome name:"+str(read_chr)
    continue


  if ( int(read_chr) < int(chr) ):
    if ( verbose):
      print "Skip reads to next chromosome:"+ str(read_chr)+ "!= "+str(chr)
    continue
      
  while ( int(read_chr) > int(chr) ):
#   print  str(read_chr)+ "!!!= "+str(chr)
   if ( not empty_bedfile):
     bedline = bedfile.readline().rstrip()
     if (not bedline):
         empty_bedfile = 1
         break

     (chr, start, end) = bedline.split("\t")

     if ( chr == "X"):
        chr = 23 

     if ( chr == "Y"):
        chr = 24 

     if ( verbose ):
        print "Next chromosome:"+read.qname+"\t"+str(read_chr)+"\t"+str(read_start)+"\t"+str( read_end )+"\t"+str(read.alen)
        print "New bed region: "+str(chr)+"\t"+str(start)+"\t"+str(end)+""
 
   if (int(read_chr) != int(chr)):
     break
  # For debugging we only look on chr 1
#  if ( int(read_chr) >= 2 ):
#    break



  # while the read start is larger than the bed block end read in new blocks.
  while ( int(read_start) > int(end) ):
   if ( not empty_bedfile):
      bedline = bedfile.readline().rstrip()
      if ( not bedline ):
          empty_bedfile = 1
          break
      (chr, start, end) = bedline.split("\t")
      if ( chr == "X"):
        chr = 23 

      if ( chr == "Y"):
        chr = 24 

      if ( verbose ):
        print "Downstream:"+read.qname+"\t"+str(read_chr)+"\t"+str(read_start)+"\t"+str( read_end )+"\t"+str(read.alen)
        print "New bed region: "+str(chr)+"\t"+str(start)+"\t"+str(end)+""

   if (int(read_chr) != int(chr)):
      break

  if ( int(read_chr) < int(chr) ):
    if (verbose):
        print "Skip reads to next chromosome:"+ str(read_chr)+ "!= "+str(chr)
    continue

  if ( (int(read_start) >= int(start) and int(read_start) <= int(end)) or
       (int(read_end)   >= int(start) and int(read_end) <= int(end))   or 
       (int(read_start) <= int(start) and int(read_end) >= int(end))      ):

    mq0 += 1

    if ( read.mapq >= 20 ):
      mq20 += 1

    if ( read.mapq >= 50 ):
      mq50 += 1


#   print read.qname

    if ( verbose ):
      print "Overlaps:"+read.qname+"\t"+str(read_chr)+"\t"+str(read_start)+"\t"+str( read_end )+"\t"+str(read.alen)
  else:
    if ( verbose ):
      print "No overlap:"+read.qname+"\t"+str(read_chr)+"\t"+str(read_start)+"\t"+str( read_end )+"\t"+str(read.alen)


print "MQ0: "+str(mq0)+", MQ20: "+str(mq20)+", MQ50: "+str(mq50)




