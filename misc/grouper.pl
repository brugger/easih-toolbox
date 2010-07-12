#!/usr/bin/perl 
# 
# groups and find the mean of depth of coverage for a region. Used when doing a chromosome plot
# 
# 
# Kim Brugger (25 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $width = shift || 10000;

my $slide = $width;

my $sum = 0;
my ($pos, $depth);
while(<>) {

  chomp;
  ($pos, $depth) = split(/\s+/, $_);

  if ( $pos > $slide) {
    print "$slide\t".($sum/$width)."\n";
    $slide += $width;
    $sum = 0;
  }
  $sum += $depth;
  
}

print "$pos\t".($sum/$width)."\n" if ( $sum && $pos);
