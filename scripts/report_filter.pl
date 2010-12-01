#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (21 Sep 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts('b:r:t:hb:dDf:F:l:', \%opts);

usage() if ( $opts{h});

my $bed_file        = $opts{b};
my $region          = $opts{r};
my $snp_file        = $opts{t} || usage();
my $in_dbsnp        = $opts{d};
my $not_dbsnp       = $opts{D};
my $filter          = $opts{g};
my $filled_field    = $opts{f};
my $unfilled_field  = $opts{F};
my $leeway          = $opts{l} || 0;

my $regions = readin_bed( $bed_file ) if ( $bed_file );

if ( $region ) {
  $region =~ s/,//g;
  my ($chr, $start, $end) = $region =~ /(.*?):(\d+)-(\d+)/;
  push @{$$regions{ $chr }}, [$start, $end];
}

open(my $snps, $snp_file ) || die "Could not open '$snp_file':$!\n";
while( <$snps> ) {

  chomp;
  s/^\s+//;

#  next if ( $in_dbsnp  && ! /\trs\d+/);
#  next if ( $not_dbsnp && /\trs\d+/);

#  next if ( $filter && !/$filter/);

  my @fields = split("\t", $_);

#  next if ( defined $filled_field  && (!$fields[$filled_field] || $fields[$filled_field] eq ""));
#  next if ( defined $unfilled_field  && $fields[$unfilled_field] &&  $fields[$unfilled_field] ne "");

  next if ( ! $fields[0]);

  my ($chr, $start, $end) = $fields[0] =~ /(.*?):(\d+)-(\d+)/;

  if ( ! $chr || !$start || ! $end ) {
    ($chr, $start) = $fields[0] =~ /(.*?):(\d+)/;
    $end = $start;
  }

  next if ( ! $chr || !$start || ! $end );

  if ( $regions ) {

    foreach my $se (@{$$regions{$chr}} ) {
      my ($Rstart, $Rend) = @$se;

      if ( ($Rstart - $leeway <= $start  && $Rend + $leeway >= $end) || 
	   ($start  - $leeway <= $Rstart && $end + $leeway >= $Rend) || 
	   ($start  - $leeway <= $Rstart && $end + $leeway >= $Rstart) ||
	   ($Rstart - $leeway <= $start  && $Rend + $leeway >= $start) ) {
	print "$_\n";
	last;
      }
      
    }
  }  
  else {
    print "$_\n";
  }
}


# 
# 
# 
# Kim Brugger (11 May 2010)
sub readin_bed {
  my ( $infile, $merge ) = @_;

  $merge = 1;
  my %res;

  open (STDIN, $infile) || die "Could not open '$infile': $!\n" if ( $infile );
  while(<STDIN>) {

    chomp;
    my ($chr, $start, $end) = split("\t", $_);

    ($chr, $start, $end) = $_ =~ /(.*?):(\d+)-(\d+)/
	if ( ! $start );
    
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
  }

  return \%res;
}



# 
# 
# 
# Kim Brugger (12 Jul 2010)
sub usage {
  
  $0 =~ s/.*\///;
  print "Filters SNPs in various ways\n";
  print "Usage: $0 -b<ed file, with regions of interest> -r[egion to extract from] -t[ab file (report, snp, indel)]\n";
  exit;

}
