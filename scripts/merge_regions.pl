#!/usr/bin/perl 
# 
# merge regions/bed entries into one.
# 
# 
# Kim Brugger (11 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my $regions = readin_bed( );

foreach my $key ( keys %$regions ) { 
  
  foreach my $region ( @{$$regions{ $key }} ) {
    
    print "$key:" . join("-", @$region) . "\n";
    
  }
}



# 
# 
# 
# Kim Brugger (11 May 2010)
sub readin_bed {
  my %res;
  my $merge = 1;
  my $flanking = 1;
  while(<STDIN>) {

    chomp;
    my ($chr, $start, $end) = split("\t", $_);

    ($chr, $start, $end) = $_ =~ /(\w+):(\d+)-(\d+)/ if ($start);
    
    
    push @{$res{$chr}}, [$start, $end];
  }

  if ( $merge ) {
    
    foreach my $key ( keys %res ) {
      
      @{$res{$key}} = sort { $$a[0] <=> $$b[0] } @{$res{$key}};
      my @tmp;
      my @data = @{$res{$key}};
      
      for(my $i=0;$i< @data; $i++) {
	
	# need at least one element in the array, so push and move on.
	if ( ! @tmp ) {
	  push @tmp, $data[ $i ];
	  next;
	}
	
	# contained in the region
	if ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ] - $flanking &&
	     $data[ $i ][ 1 ] <= $tmp[ -1 ][ 1 ] + $flanking) {
	  next;
	}
	# overlapping
	elsif ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ] - $flanking  &&
		$data[ $i ][ 0 ] <= $tmp[ -1 ][ 1 ] + $flanking) {
	  
	  $tmp[ -1 ][ 1 ] = $data[ $i ][ 1 ];
	}
	# There is a gap between the end block and this one. Just push it on the end of  the array!
	else {
	  push @tmp, $data[ $i ];
	}
      }
      @{$res{$key}} = @tmp;
    }
  }

  return \%res;
}
