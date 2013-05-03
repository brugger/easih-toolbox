#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (15 Oct 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my $opts = 'r:';
my %opts;
getopts($opts, \%opts);


my $infile1 = shift;
my $infile2 = shift;

my $reference = "/scratch/kb468/GATK_bundle/b37/human_g1k_v37.fasta";

$reference = $opts{r} if ($opts{r});

my $cmd  = "bwa_0.6.1-tpx  aln -t 8 -e 5 -q 15 -f $infile1.sai $reference $infile1";
print "$cmd\n";
system $cmd;

$cmd  = "bwa_0.6.1-tpx  aln -t 8 -e 5 -q 15 -f $infile2.sai $reference $infile2";
system $cmd;


$cmd  = "bwa_0.6.1-tpx sampe -t 8 $reference $infile1.sai $infile2.sai $infile1 $infile2 | samtools view -Sb - | samtools sort - $infile1";
print "$cmd\n";
system $cmd;

$cmd  = "samtools index $infile1.bam";
system $cmd;

$cmd  = "samtools flagstat $infile1.bam > $infile1.flagstat";
system $cmd;


exit;

__END__

$cmd  = "/scratch/kb468/pysam-0.6/Exome_OT.py -b $infile1.bam -B /data/B10/capture.bed > $infile1.OT";
system $cmd;

$cmd  = "/scratch/kb468/pysam-0.6/Exome_OT.py -b $infile.bam -B /data/B10/capture.bed > $infile.OT";
system $cmd;



