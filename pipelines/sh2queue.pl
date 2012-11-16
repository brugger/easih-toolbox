#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (10 Oct 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

while(<>) {
  chomp;
  system "echo \" $_ \" | qsub -cwd\n"
}
