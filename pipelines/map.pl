#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (15 Oct 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $infile = shift;

my $reference = "/scratch/kb468/GATK_bundle/b37/human_g1k_v37.fasta";

my $cmd  = "bwa_0.6.2-tpx  aln  -e 50 -f $infile.sai $reference $infile";
system $cmd;


$cmd  = "bwa_0.6.2-tpx samse   $reference $infile.sai $infile | samtools view -Sb - | samtools sort - $infile";
system $cmd;

$cmd  = "samtools index $infile.bam";
system $cmd;

$cmd  = "samtools flagstat $infile.bam > $infile.flagstat";
system $cmd;



$cmd  = "/scratch/kb468/pysam-0.6/Exome_OT.py -b $infile.bam -B /data/B10/capture.bed > $infile.OT";
system $cmd;

$cmd  = "/scratch/kb468/pysam-0.6/Exome_OT.py -b $infile.bam -B /data/B10/capture.bed > $infile.OT";
system $cmd;



