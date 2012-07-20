#!/usr/bin/perl 
# 
# pipeline for few genes analysis, mainly for clinical use.
# 
# 
# Kim Brugger (17 Jan 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

BEGIN {
  use vars qw/$path/; 
  $path = $0;
  if ($path =~ /.*\//) {
    $path =~ s/(.*\/).*/$1/;
  }
  else {
    $path = "./";
  }
  push @INC, $path;
}

use Getopt::Std;

#use lib '/home/cjp64/git/easih-pipeline/modules';
use lib '/home/kb468/easih-pipeline/modules';
use lib '/home/kb468/easih-toolbox/modules';


use EASIH::Misc;

while(<>) {
  if (/^\@/) {
    s/(.*SN:\d*):\d+-\d+(.*)/$1$2/ if ( /SN:\d+/);
    print  $_ if ( ! /SQ/);
  }
  else {
    my @F = split("\t");
    if ( $F[2] =~ /(.*?):(\d+)-(\d+)/) {
      my ($chr, $start, $end) = ($1,$2,$3);
      $F[3] += $start - 1;
      $F[2] = $chr;
      print join("\t", @F);
    }
  }
}




