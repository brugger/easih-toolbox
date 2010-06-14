#!/usr/bin/perl 
# 
# On-Off target scores
# 
# 
# Kim Brugger (11 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


my %lengths;

my ($on_target, $off_target) = (0,0);

my %res = ();

while(<>) {
  
  chomp;

  my ($id, $flag, $chr, $start, undef, $cigar, undef, undef, undef, $read, undef) = split("\t");
  

  my (@parts) =  $cigar =~ /(\d*\w)/g;

  my $map_length = 0;
  
  foreach my $part ( @parts ) {
    
    my ($length, $type) =  $part =~ /(\d*)(\w)/;
    $length ||= 1;
    
    next if ( $type ne "M");
    
    $map_length += $length;
  }


  $res{ $chr }{ $map_length }++;
		    
}

foreach my $chr ( keys %res ) {
  
  my ($count, $sum) = (0,0);

  foreach my $length ( keys  %{$res{ $chr }} ) {

    $count += $res{ $chr }{ $length };
    $sum   += $length * $res{ $chr }{ $length };
  }

  printf("$chr: %5.2f \n", $sum/$count);
}
