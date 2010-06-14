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

my $bam_file    = shift || die "No bamfile given!\n";
my $no_cov_file = shift || undef;
my %exon_res;

my @no_cov;

while(<>) {


  $_ =~ s/\r//g;
  $_ =~ s/\n//g;
  my ($chr, $start, $end) = split("\t", $_);
  
  $chr = "chr$chr" if ( $chr !~ /chr/);
  
  my $region = "$chr:$start-$end";

  my ($summed_depth, $length) = (0,0);
  
  open (my $bam_pipeline, " $bam2depth $bam_file $region | ") || die "Could not open bam2depth pipeline: $! ($bam2depth $bam_file $region)\n";

  my $last_pos;
  while ( <$bam_pipeline> ) {
    chomp;
    my( $region, $pos, $level) = split("\t");

    $summed_depth += $level;

    if ( $last_pos &&  $last_pos + 1 < $pos ) {
      push @no_cov, "$region:" .($last_pos + 1 )."-" .($pos -1 )."\n";
    }

    $last_pos = $pos;
  }

  if (! $last_pos ||  $last_pos +1 < $end ) {
    $last_pos ||= $start;
    push @no_cov, "$region:" .($last_pos + 1 )."-" .($end )."\n";
  }


  my $mean = int($summed_depth/ ($end- $start + 1));

  print  "$region\t$mean\n";
}



if ( $no_cov_file && @no_cov ) {
  
  open (my $out, "> $no_cov_file") || die "Could not open '$no_cov_file':$!\n";
  print $out "@no_cov";
  close ($out);

}
