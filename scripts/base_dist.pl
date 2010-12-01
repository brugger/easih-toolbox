#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (28 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my $bam = shift;


while(my $region = shift || <>) {
  chomp($region);

  if ( $region =~ /(\w+):(\d+)-(\d+)/ || $region =~ /(\w+)\t(\d+)\t(\d+)/) {
    my ($chr, $start, $end) = ($1, $2, $3);
    
    for (my $i = $start; $i <= $end; $i++) {
      base_dist( $1, $i );
    }
  }
  elsif ($region =~ /(\w+):(\d+)/ || $region =~ /(\w+)\t(\d+)/) {
    my ($chr, $pos) = ($1, $2);
    base_dist( $chr, $pos);
  }
  else {
    print STDERR "Does not know how to parse region and position from $region\n";
  }
}


# 
# 
# 
# Kim Brugger (29 Apr 2010)
sub base_dist {
  my ( $chr, $SNP_pos) = @_;

  if ( ! $bam ) {
    print STDERR "need a bam file for finding base distribution\n";
    exit;
  }

  $chr =~ s/chr//;
#  $chr = "chr$chr";

  my %base_stats = ( A => 0, C => 0, G => 0, T => 0, N => 0);
  my %pos_stats;
  my %qual_stats;
  my $total = 0;

  open (my $st_pipe, "samtools view -q 5 $bam $chr:$SNP_pos-$SNP_pos | ") || die "Could not open samtools pipe: $!";

#  print "samtools view -q 5 $bam $chr:$SNP_pos-$SNP_pos\n";

  while(<$st_pipe>) {
    chomp;
    my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");

    ($sequence, $quality) = patch_alignment($sequence, $quality, $cigar);

    my @seq  = split("", $sequence);
    my @qual = split("", $quality);    

    my $base_pos = $SNP_pos - $pos;
    my $base = $seq[ $base_pos ];

    if ( !$base ) {
      print "$SNP_pos $_\n \n";
      print STDERR "FAILED !!!! \n";
      exit;
    }

    my $qual = ord($qual[ $base_pos ])-33;
    $base_stats{ $base }++;
    $pos_stats{$base}{$pos}++;
    push @{$qual_stats{$base}},  $qual;
    $total++;
  }


  if (0) {
    foreach my $base ('A', 'C', 'G', 'T' ) {
      my $dups = keys %{$pos_stats{ $base}};
      if ( $dups == 1) {
	$total -= $base_stats{ $base };
	$base_stats{$base} = 0;
	delete $qual_stats{$base};
      }
    }    
  }

  foreach my $key ( keys %qual_stats ) {
    my $sum = eval join '+', @{$qual_stats{$key}};
    my $count = @{$qual_stats{$key}};
    my $avg_qual      =  sprintf("%.2f",$sum / $count);
    $avg_qual ||= 0;
    $qual_stats{$key} = $avg_qual;
  }

  print "$chr:$SNP_pos\t";
  if ( $total ) {
    my %res;
#  foreach my $base (sort {$base_stats{$b} <=> $base_stats{$a}} keys %base_stats ) {
    foreach my $base ('A', 'C', 'G', 'T' ) {
      my $perc = sprintf("%.2f", $base_stats{$base}/$total*100);
      my $qual = $qual_stats{$base} || 0;
      my $dups = keys %{$pos_stats{ $base}};
      print "$base: $base_stats{$base}($perc%)/$qual/$dups\t";
    }
  }
  print "\n";
#  exit;
#  print "done\n";
}


# 
# 
# 
# Kim Brugger (20 Jul 2009)
sub patch_alignment {
  my ( $seq, $qual, $cigar ) = @_;

  return ($seq, $qual) if ( $cigar !~ /[DIS]/);
  
  my @seq  = split("", $seq );
  my @qual = split("", $qual );


  my (@cigar) = $cigar =~ /(\d*\w)/g;

  my $offset = 0;

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)



  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d+)(\w)/;

    if ( $type eq 'M') {
      $offset += $length;
      next;
    }
    elsif ( $type eq "D") {
      my @dashes = split("", "-"x$length);
      splice(@seq,  $offset, 0, @dashes);
      splice(@qual, $offset, 0, @dashes);
      $offset += $length;
    }
    elsif ( $type eq "I" || $type eq "S" ) {
      splice(@seq,  $offset, $length);
      splice(@qual, $offset, $length);
    }    

  }


  return (join("", @seq), join("",@qual));
}
