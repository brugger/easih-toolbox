#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (24 Jun 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


use lib '/home/kb468/projects/easih-flow/modules';
use lib '/home/kb468/easih-pipeline/modules';
use EASIH::JMS;

my $freeze_file = shift || usage();


EASIH::JMS::restore_state($freeze_file);
print EASIH::JMS::report();    
print EASIH::JMS::total_runtime() . "\n";

my $solid2fq  = '/home/kb468/bin/solid2fastq.pl';
my $bwa       = '/home/easih/bin/bwa';
my $samtools  = '/home/easih/bin/samtools';
my $fq_split  = '/home/kb468/bin/fastq_split.pl';

my $extra_report = "";
$extra_report .= " Binaries used..\n";
$extra_report .= `ls -l $samtools`;
$extra_report .= `ls -l $bwa` . "\n";

print "$extra_report";

print EASIH::JMS::full_report();    




# 
# 
# 
# Kim Brugger (24 Jun 2010)
sub usage {
  
  $0 =~ s/.*\///;
  die "USAGE: $0 pipeline-file\n";
  
}

