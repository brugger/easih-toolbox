#!/usr/bin/perl 
# 
# On-Off target scores
# 
# 
# Kim Brugger (11 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;


my %opts;
getopts('b:B:f:h', \%opts);
usage() if ( $opts{h});

my $bed_file = $opts{B} || shift;
my $bam_file = $opts{b} || usage();
my $flanking = $opts{f} || 100;

my $regions = readin_bed( $bed_file, 1 );

my $samtools = `which samtools`;
chomp( $samtools);

my $all_on_target = 0;

foreach my $chr ( keys %$regions ) {

  my $on_target = 0;
  
  foreach my $se ( @{$$regions{$chr}}) {
    
    $chr =~ s/chr//i;
#    $chr ="chr$chr";

    my ($start, $end) = @$se;
    $start = $start - $flanking;
    $end   = $end + $flanking;

    my $region = "$chr:$start-$end";
    open (my $st_pipe, "$samtools view -F 0x0404 $bam_file $region | ") || die "Could not open samtools pipe: $!\n";

    while(<$st_pipe>) {
      $on_target++;
    }
  }

  $all_on_target += $on_target;

  print "$chr on target: $on_target\n";
}




# 
# 
# 
# Kim Brugger (11 May 2010)
sub readin_bed {
  my ( $infile, $merge ) = @_;

  my %res;

  open (STDIN, $infile) || die "Could not open '$infile': $!\n" if ( $infile );
  while(<STDIN>) {

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



# 
# 
# 
# Kim Brugger (12 Jul 2010)
sub usage {
  
  $0 =~ s/.*\///;
  print "Finds the nr of reads that can be put onto a target (on a chromosome basis). Overlapping regions are merged.\n";
  print "Usage: $0 -b<am file> -f[lank, default 100] -B[ed file]/bedfile/STDIN\n";
  exit;

}
