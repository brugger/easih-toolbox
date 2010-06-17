#!/usr/bin/perl 
# 
# Converts a fastq file to a fasta file.
# 
# 
# Kim Brugger (24 Mar 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $fastq_file = "";
my $fasta_file = "";

&GetOptions(
            'fastq_file'  => \$fastq_file,
            'fasta_file'  => \$fasta_file,
           );

$fastq_file ||= shift;

if (! $fastq_file ) {
  system "perldoc $0";
  exit;
}

my $basename = basename($fastq_file, qw/.fastq .fq/);
$fasta_file = $basename . ".fa";


open ( my $fastq, "$fastq_file" ) || die "Could not open '$fastq_file': $!\n";
while(my $header = <$fastq> ) {
  my $sequence   = <$fastq>;
  my $strand     = <$fastq>;
  my $quality    = <$fastq>;

  $header =~ s/^\@/\>/;

  print $header. nice_fasta($sequence);
#  exit;
}


sub nice_fasta {
 my ($seq) =  @_;
 my $j = 0;
 my $count = length $seq;                                                       
 $seq =~ s/\n//g;
 my $line_length = 60;
 my $res = "";

 while ($j < $count) {
   $res .= substr($seq, $j, $line_length). "\n";
    $j += $line_length;
 } 

 return $res;
}



=pod

=head1 NAME

fastq_trim.pl trims a fastq file based on quality and/or sequence length

=head1 OPTIONS

Usage: fastq_trim.pl <options>.

=over

=item B<-fastq_file F<fastq_file>>: 
    
The fastq file to analysse


=item B<score F<cutoff score>>

Min score to filter on.


=item B<length F<min length>>

sequences shorter than this length is dropped. If the sequence is trimmed as 
well this step is done after the trimming.


=cut
