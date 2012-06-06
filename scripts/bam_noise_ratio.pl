#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (06 Jun 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $bamfile    = shift;
my $reference = shift || "";

Usage() if ( ! $bamfile );


my ($mapped, $deletions, $insertions) = (0,0,0);
my ($mapped_block, $deletion_block, $insertion_block) = (0,0,0);


my $max_reads = -1;
$max_reads = 20000;

open(my $bam_ph, "samtools view $bamfile | ") || die "Could not open bam-stream 'samtools view $bamfile': $!\n";
while (<$bam_ph>) {
  chomp;

  my @F = split("\t");
  # unmapped read...
  next if ($F[ 1 ]  & 0x0004 );
  my $cigar = $F[ 5 ];

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)

  foreach my $part ( $cigar =~ /(\d*\w)/g) {
    my ($length, $type) =  $part =~ /(\d*)(\w)/;
    $length ||= 1;

    next if ( $type eq "S" || $type eq "N" || $type eq "H" || $type eq "P");
    $mapped     += $length if ( $type eq "M" );
    $deletions  += $length if ( $type eq "D" );
    $insertions += $length if ( $type eq "I" );

    $mapped_block++     if ( $type eq "M" );
    $deletion_block++   if ( $type eq "D" );
    $insertion_block++  if ( $type eq "I" );

  }
#  printf("%.2f $max_reads\n", ($deletions+$insertions)/$mapped);
  last if ( --$max_reads == 0 );

}		      

#printf("Indel ratio of mapped bases: %.2f%%\n", ($deletions+$insertions)*100/$mapped);
#printf("Indel frequency: %.2f%%\n", ($deletion_block+$insertion_block)/($mapped_block));

printf("Noise Index: %.2f\n", (($deletions+$insertions)*100/$mapped)*(($deletion_block+$insertion_block)/$mapped_block));




# 
# 
# 
# Kim Brugger (06 Jun 2012)
sub Usage {
  $0 =~ s/.*\///;

  print "USAGE: Calculates a noise index for alignened bases in a bam file. The index is defined as the ratio indel vs aligned bases times the frequency of indel blocks\n";
  print "USAGE: $0 <BAMfile>\n" ;
      
  exit -1;
}

  



