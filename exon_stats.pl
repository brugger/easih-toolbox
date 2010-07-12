#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (25 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use Getopt::Std;


my %opts;
getopts('b:B:n:h', \%opts);
usage() if ( $opts{h});

my $bam2depth = `which bam2depth` || "/usr/local/bin/bam2depth";
chomp($bam2depth);

my $bed_file    = $opts{B} || shift;
my $bam_file    = $opts{b} || usage();
my $no_cov_file = $opts{n} || undef;
my %exon_res;

my @no_cov;

my $regions = readin_bed( $bed_file);

foreach my $chr ( keys %$regions ) {
  foreach my $se ( @{$$regions{$chr}}) {

    my ($start, $end) = @$se;
    $chr = "chr$chr" if ( $chr !~ /chr/);

    my $st_region = "$chr:$start-$end";

    my ($summed_depth, $length) = (0,0);
  
    open (my $bam_pipeline, " $bam2depth $bam_file $st_region | ") || die "Could not open bam2depth pipeline: $! ($bam2depth $bam_file $st_region)\n";

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




# 
# 
# 
# Kim Brugger (11 May 2010)
sub readin_bed {
  my ( $infile ) = @_;

  open (STDIN, $infile) || die "Could not open '$infile': $!\n" if ( $infile );

  my %res;

  while(<STDIN>) {

    $_ =~ s/\r//g;
    $_ =~ s/\n//g;
    my ($chr, $start, $end) = split(/\s+/, $_);
    
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

  return \%res;
}




# 
# 
# 
# Kim Brugger (12 Jul 2010)
sub usage {
  
  $0 =~ s/.*\///;
  print "Finds the mean depth of exons/targets plus can report regions with no/partial coverage.\n";
  print "Usage: $0 -b<am file> -n[o coverage (outfile)] -B[ed file]/bedfile/STDIN\n";
  exit;

}
