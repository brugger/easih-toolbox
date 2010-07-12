#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (25 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $bam2depth = "/usr/local/bin/bam2depth";
   $bam2depth = "/home/easih/bin/bam2depth";

my $bam_file    = shift || die "No bamfile given!\n";
my $no_cov_file = shift || undef;
my %exon_res;

my @no_cov;

my %res;
while(<>) {

  chomp;
  s/\n//g;
  s/\r//g;
  my ($chr, $start, $end) = split("\t", $_);
  
  push @{$res{$chr}}, [$start, $end];
}


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
    if ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ] &&
	 $data[ $i ][ 1 ] <= $tmp[ -1 ][ 1 ]) {
      next;
    }
    # overlapping
    elsif ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ]  &&
	    $data[ $i ][ 0 ] <= $tmp[ -1 ][ 1 ]) {
      
      $tmp[ -1 ][ 1 ] = $data[ $i ][ 1 ];
    }
    # There is a gap between the end block and this one. Just push it on the end of  the array!
    else {
      push @tmp, $data[ $i ];
    }
  }
  @{$res{$key}} = @tmp;
}



foreach my $chr ( keys %res ) {

  
  foreach my $region ( @{$res{$chr}}) {

    $chr =~ s/chr//i;
    $chr ="chr$chr";

    my $start = $$region[0];
    my $end   = $$region[1];

    my $st_region = "$chr:$start-$end";

    my ($summed_depth, $length) = (0,0);
  
    open (my $bam_pipeline, " $bam2depth $bam_file $st_region | ") || die "Could not open bam2depth pipeline: $! ($bam2depth $bam_file $region)\n";

    my $last_pos;
    while ( <$bam_pipeline> ) {
    chomp;
    my( $region, $pos, $level) = split("\t");

    $summed_depth += $level;
    
    if ( $last_pos &&  $last_pos + 1 < $pos ) {
      push @no_cov, "$st_region:" .($last_pos + 1 )."-" .($pos -1 )."\n";
    }
    
    $last_pos = $pos;
  }
    
    if (! $last_pos ||  $last_pos +1 < $end ) {
      $last_pos ||= $start;
      push @no_cov, "$st_region:" .($last_pos + 1 )."-" .($end )."\n";
    }
    
    
    my $mean = int($summed_depth/ ($end- $start + 1));
    
    print  "$st_region\t$mean\n";
  }
}


if ( $no_cov_file && @no_cov ) {
  
  open (my $out, "> $no_cov_file") || die "Could not open '$no_cov_file':$!\n";
  print $out "@no_cov";
  close ($out);

}
