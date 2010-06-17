#!/usr/local/bin/perl

# Produce a fastq file from a SAM file

use strict;
use warnings;

my $VERSION = '20091221.01';

sub process;

process;

exit;

# ----------------------------------------------------------------------
sub process {

  while (my $line = <>) {    # sets $line

    chomp $line;

    next if $line =~ /^\@/xms;           # skip header lines

    my ($name, $flag, $read, $qual) = (split /\t/, $line)[0,1,9,10];

    if ($flag & 0x0010) {                # SAM revcomps reads aligned to reverse strand
      $read = reverse $read;
      $read =~ tr/ACGTacgt/TGCAtgca/;
      $qual = reverse $qual;
    }

    print "\@$name\n$read\n+\n$qual\n";

  }

  return;

}

