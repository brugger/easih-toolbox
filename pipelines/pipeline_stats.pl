#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (24 Jun 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


use lib '/home/kb468/scratch/BRC/easih-pipeline/modules';
use EASIH::Pipeline;

my $freeze_file = shift || usage();


EASIH::Pipeline::restore_state($freeze_file);
print EASIH::Pipeline::report();    
print EASIH::Pipeline::total_runtime();
print EASIH::Pipeline::real_runtime(). "\n";

my $solid2fq  = '/home/kb468/bin/solid2fastq.pl';
my $bwa       = '/home/easih/bin/bwa';
my $samtools  = '/home/easih/bin/samtools';
my $fq_split  = '/home/kb468/bin/fastq_split.pl';

print EASIH::Pipeline::full_report();    




# 
# 
# 
# Kim Brugger (24 Jun 2010)
sub usage {
  
  $0 =~ s/.*\///;
  die "USAGE: $0 pipeline-file\n";
  
}

