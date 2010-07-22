#!/usr/bin/perl 
# 
# Makes fastq files from 454 fna and qual files.
# 
# 
# Kim Brugger (24 Mar 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use File::Basename;

my $fasta_file = shift;
my $qual_file  = shift;
my $basename = basename($fasta_file, qw/.fasta .fna/);
my $fastq_file = $basename . ".fastq";

$qual_file ||= $basename . ".qual";

my %seqs;

$/ = ">";

open (my $fasta, " $fasta_file") || die "Could not open '$fasta_file': $!\n";
while (<$fasta>) {
  $_ =~ s/\r//g;
  chomp;
  my ($header, @sequence) = split /\n/;
  next if ( ! $header );
  $header =~ s/^(\w+) .*/$1/;
  $seqs{ $header } = join('', @sequence);
}
close( $fasta);

open (my $qual, "$qual_file") || die "Could not open '$qual_file': $!\n";
open (my $fastq, "> $fastq_file") || die "Could not open '$fastq_file': $!\n";


while (<$qual>) {
  $_ =~ s/\r//g;
  chomp;
  my ($header, @qualities) = split /\n/;
  next if ( ! $header );
  @qualities = split(" ",  join ' ', @qualities);
  print $fastq "\@"."$header\n$seqs{$header}\n+\n";
  my $fastq_line = "";
  map { $fastq_line .= chr($_ + 33) } @qualities;
  print $fastq "$fastq_line\n";
}

close( $qual );
close( $fastq );
