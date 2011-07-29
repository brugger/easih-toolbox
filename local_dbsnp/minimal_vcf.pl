#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (17 Jan 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

while(<>) {
  if (/#/) {
    print;
    next;
  }
  else {
    chomp;
    my @f = split("\t");
    @f = @f[0..6];
    print join("\t", @f) . "\n";
  }
}
