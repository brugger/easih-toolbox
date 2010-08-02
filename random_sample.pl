#!/usr/bin/perl 
# 
# Testing code for randomly sampling XXX MB of data from a fq file.
# 
# 
# Kim Brugger (02 Aug 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my %used;

my $infile = shift;


my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($infile);

open (my $file, $infile) || die "Could not open '$infile': $!\n";
my $read = 0;
my $goal = 1*1024;#000;
while( $goal > $read ) {

  my $random_pos = int(rand($size));
  print "going to $random_pos\n";
  seek( $file, $random_pos, 0);
  while ( <$file> ) {
    last if ( $_ =~ /^\@/);
  }

  if ($_ &&  $_ =~ /^\@/ && !$used{ $_ }) {
    my $name = $_;
    my $seq  = <$file>;
    my $str  = <$file>;
    my $qual = <$file>;

    print "$name\n";
    $read += length("$name$seq$str$qual");
    $used{$file}++;
  }


}

