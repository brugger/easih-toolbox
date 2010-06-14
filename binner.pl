#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (25 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


my %res;

my $width = shift || 10;

while(<>) {

  chomp;
  my ($value, $depth) = split(/\s+/, $_);
  
  my $bin_value = int($value/$width)*$width;
  
  $res{ $bin_value } += $depth;
  
}


foreach my $key ( sort {$a <=>$b} keys %res ) {
  print "$key\t$res{$key}\n";
}
