#!/usr/bin/perl 
# 
# On-Off target scores
# 
# 
# Kim Brugger (11 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $regions_file  = shift;
my $bam_file      = shift;
my $leeway        = shift || 100;

my $regions = readin_bed( $regions_file, 1 );

my $samtools = `which samtools`;
chomp( $samtools);

foreach my $chr ( keys %$regions ) {

  my $on_target = 0;
  
  foreach my $region ( @{$$regions{$chr}}) {

    my $start = $$region[0] - $leeway;
    my $end   = $$region[1] + $leeway;

    my $st_region = "$chr:$start-$end";
    open (my $st_pipe, "$samtools view $bam_file $st_region | ") || die "Could not open samtools pipe: $!\n";

    while(<$st_pipe>) {
      $on_target++;
    }
  }

  print "$chr on target: $on_target\n";
}




# 
# 
# 
# Kim Brugger (11 May 2010)
sub readin_bed {
  my ( $infile, $merge ) = @_;

  my %res;

  open ( my $in, $infile) || die "Could not open '$infile': $!\n";
  while(<$in>) {

    chomp;
    my ($chr, $start, $end) = split("\t", $_);
    
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
	if ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ] - $leeway &&
	     $data[ $i ][ 1 ] <= $tmp[ -1 ][ 1 ] + $leeway) {
	  next;
	}
	# overlapping
	elsif ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ] - $leeway  &&
		$data[ $i ][ 0 ] <= $tmp[ -1 ][ 1 ] + $leeway) {
	  
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
