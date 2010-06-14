#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (23 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/home/kb468/easih-misc/modules/';
use EASIH::Profile;

my $infile  = shift || "/home/kb468/X_probes.bed";
my $bamfile = shift || die "no bam file specified\n";

my %doc = ();

my $counter = 0;
open (my $in, $infile) || die "Could not open '$infile': $!\n";
while (<$in>) {
  chomp;
  
  my ($chr, $start, $end) = split("\t", $_);
  my $region = "$chr:$start-$end";

  my $profile = EASIH::Profile->New;

  $profile->from_bam($bamfile, "$region");
  
#  $profile->dump_profile;

  my $mean_coverage = $profile->mean_coverage();
#  print "MC : $mean_coverage \n";
  $doc{ int($mean_coverage)}++;
#  last if ( ++$counter >= 100);
}


foreach my $key ( sort { $a <=> $b }keys %doc ) {
  print "$key\t$doc{ $key}\n";
}
