#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (23 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

use lib '/home/kb468/easih-misc/modules/';
use EASIH::Profile;


my %opts;
getopts('b:B:h', \%opts);
usage() if ( $opts{h});


my $bed_file = $opts{B} || shift;
my $bam_file = $opts{b} || usage();
my $flanking = $opts{f} || 200;

my %doc = ();

my $counter = 0;
open (my $in, $bed_file) || die "Could not open '$bed_file': $!\n";
while (<$in>) {
  chomp;
  
  my ($chr, $start, $end) = split("\t", $_);
  my $region = "$chr:$start-$end";

  my $profile = EASIH::Profile->New;

  $profile->from_bam($bam_file, "$region");
  

  my $mean_coverage = $profile->mean_coverage();
  $doc{ int($mean_coverage)}++;
}


foreach my $key ( sort { $a <=> $b }keys %doc ) {
  print "$key\t$doc{ $key}\n";
}



# 
# 
# 
# Kim Brugger (12 Jul 2010)
sub usage {
  
  $0 =~ s/.*\///;
  print "Finds mean probe coverage, output is tab seperated mean_depth and nr of regions with that depth\n";
  print "Usage: $0 -b<am file> -B[ed file]\n";
  exit;

}
