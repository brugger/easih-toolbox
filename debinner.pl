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
  my ($value, $count) = split(/\s+/, $_);

  for(my $i = 0; $i < $count; $i++) {
    print "$value\n";
  }
}

