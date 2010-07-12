#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Mar 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use POSIX qw(ceil floor);
use Getopt::Std;


my %opts;
getopts('b:B:f:h', \%opts);
usage() if ( $opts{h});

my $bam2depth = `which bam2depth` || "/usr/local/bin/bam2depth";
chomp($bam2depth);

my $bed_file = $opts{B} || shift;
my $bam_file = $opts{b} || usage();
my $flanking = $opts{f} || 200;



my $regions = readin_bed( $bed_file);


foreach my $chr ( keys %$regions ) {
  foreach my $se ( @{$$regions{$chr}}) {
    
    my ($start, $end) = @$se;
    $chr = "chr$chr" if ( $chr !~ /chr/);
  
    stats($chr, $start - $flanking, $start, $start);
    stats($chr, $end, $end + $flanking, $end);
  }
}


# 
# 
# 
# Kim Brugger (14 Jun 2010)
sub stats {
  my ( $chr, $start, $end, $count_start ) = @_;  
  
  my $region = "$chr:$start-$end";

  my ($summed_depth, $length) = (0,0);
  
  open (my $bam_pipeline, " $bam2depth $bam_file $region | ") || die "Could not open bam2depth pipeline: $! ($bam2depth $bam_file $region)\n";
  while ( <$bam_pipeline> ) {
    chomp;
    my( $region, $pos, $level) = split("\t");

    next if ( $pos == $end || $pos == $start );

    if ( $region ) {
      print "". ($pos - $count_start )."\t$level\n";
    }
  }
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
  print "Finds the depth surrounding a target.\n";
  print "Usage: $0 -b<am file> -f[lank, default 200] -B[ed file]/bedfile/STDIN\n";
  exit;

}

