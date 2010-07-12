#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Mar 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
no warnings 'recursion';
use Data::Dumper;
use POSIX qw(ceil floor);

my $bam2depth = "/usr/local/bin/bam2depth";

my $bam_file = shift || die "No bamfile supplied\n";
my $flanking = shift || 200;

while (<>) {

  $_ =~ s/\r//g;
  $_ =~ s/\n//g;
  my ($chr, $start, $end) = split("\t", $_);
  
  $chr = "chr$chr" if ( $chr !~ /chr/);
  
  stats($chr, $start - $flanking, $start, $start);
  stats($chr, $end, $end + $flanking, $end);


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

