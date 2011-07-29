#!/usr/bin/perl -w
#*************************************************************************
#
#   File:       KickoffAnalysisOnHPC.pl
#   
#   Version:    V1.0
#   Date:       31.05.11
#   Function:   Uploading sequencing data to HPC
#   
#   Copyright:  (c) University of Cambridge / S.V.V. Deevi 2011
#   Author:     S.V.V. Deevi
#   Address:    Eastern Sequence Informatics Hub (EASIH),
#               The University of Cambridge,
#               Level 4, Laboratory Block Box 268,
#               Addenbrooke's Hospital,
#               Cambridge,
#               CB2 0QQ,
#               UK.
#   EMail:      svvd2@cam.ac.uk (or) sri.deevi@easih.ac.uk
#               
#
#   Description:
#   ============
#   Script for kicking off analysis (MARIS) on HPC. 
#
#*************************************************************************
#
#   Usage:
#   ====== 
#   $0 -h<help> -p<platform:solid or illumia> -r<reference(default:GRCHG37)> -d<dbsnp(default:dbsnp_130_b37.rod)> [input directory(ideal: /scratch/<project name>), default = current directory]
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  31.05.11 Original
#*************************************************************************

use strict;
use Getopt::Std;

my %opts;
getopts('hp:r:d:', \%opts);
my $platform          = $opts{p} || usage();
my $reference         = $opts{r} || 0;
my $dbsnp             = $opts{d} || 0;

my $indir = shift;

usage(), if($opts{h});

$platform = uc($platform), if($platform);

$dbsnp = "/scratch2/easih/refs/GATK/dbsnp_130_b37.rod", if(!($dbsnp));    

if($platform eq "SOLID")
{
    $reference = "/scratch2/easih/refs/human_1kg/bwa_cs/human_g1k_v37.fasta", if(!($reference));    
}
elsif($platform eq "ILLUMINA")
{
    $reference = "/scratch2/easih/refs/human_1kg/bwa/human_g1k_v37.fasta", if(!($reference));    
}


my $pwd;

if($indir)
{
    chomp $indir;
    $indir =~ s/\/$//;
    $pwd = "$indir";
}

chomp($pwd = `pwd`), unless($pwd);
print "\n Input directory (default = current directory): $pwd \n";

my(@fqfiles) = PickRecentfqFile($pwd);


foreach(@fqfiles)
{
    if(/(.*\/[A-Z]\d{7}.*).1.fq.*/)
    {
	my $basename = $1;
	
	system("nohup /home/kb468/easih-toolbox/pipelines/MARIS.pl -Q $basename -R $reference -d $dbsnp -p $platform &");
    }
}    


#################################################################################
sub PickRecentfqFile
{
    my($pwd) = @_;
    
    chomp(my $test =  `ls -rt ${pwd}/*fq.gz | tail -1`); 
    print "\n Last updated file in the directory: $test \n";

    chomp(my $recent = `ls -lrt --time-style=long-iso  ${pwd}/*fq.gz | awk '{print \$6}' | tail -1`); 
    print "\n date of recent fastq file (fastq files which share same date will be pulled out for further processing): $recent \n";

    chomp(my $today = `date '+%F'`);
    
    chomp(my $date_diff = `echo "\$((\$(date -d '$today' "+%s") - \$(date -d '$recent' "+%s")))" | awk '{printf "%.0f", \$1/(24*60*60)}'`);
    
    print "\n\n *** ALERT! Files of interest are older than a day: Date difference is $date_diff day(s) ***\n\n", if(($date_diff>1));

    chomp(my $list = `ls -lrt --time-style=long-iso ${pwd}/*fq.gz | grep "$recent" | awk '{print \$8}' | sort | xargs  echo`);
    print "\n List of fastq files which share same date as the recent file in the directory: $list \n";
    
    my(@files) = split(" ",$list);

    return(@files);
}


#################################################################################
sub usage
{
  $0 =~ s/.*\///;
  print STDERR "\n\nDescription: Kicksoff analysis (MARIS) on HPC\n";
  print STDERR "\n\nDetails: Picks up recent fastq files in the directory and submits analysis(MARIS) jobs on HPC\n";
  print STDERR "\nUsage: $0 -h<help> -p<platform:solid or illumia> -r<reference(default:GRCHG37)> -d<dbsnp(default:dbsnp_130_b37.rod)> [input directory(ideal: /scratch/<project name>), default = current directory]\n\n\n";
  exit;
}
