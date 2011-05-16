#!/usr/bin/perl 
# 
# program to trim off low quality sequence from the ends of 
# (454) sequences
# 
# Kim Brugger (24 Mar 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;



my $fastq_file = "";
my $min_score  = 0;
my $length     = 0;


&GetOptions(
            'fastq_file'  => \$fastq_file,
            'score:n'     => \$min_score,
            'length:n'    => \$length,
           );

$fastq_file ||= shift;
$min_score  += 33;

if (! $fastq_file ) {
  system "perldoc $0";
  exit;
}

open ( my $fastq, "$fastq_file" ) || die "Could not open '$fastq_file': $!\n";
while(my $header = <$fastq> ) {
  my $sequence   = <$fastq>;
  my $strand     = <$fastq>;
  my $quality    = <$fastq>;

  my @seqs  = split("", $sequence );
  my @quals = split("", $quality );
  

  if ( $min_score > 33 ) {

    my $start_trim = 0;
    for (; $start_trim < @quals; $start_trim++) {
      last if ( ord $quals[ $start_trim ] > $min_score);
    }

    my $end_trim = @quals - 1;
    for (; $end_trim > $start_trim; $end_trim--) {
      last if ( ord ($quals[ $end_trim ]) > $min_score);
    }
    
    $sequence = substr( $sequence, $start_trim, $end_trim - $start_trim + 1);
    $quality  = substr( $quality , $start_trim, $end_trim - $start_trim + 1);
  }

  next if ( $length && $length > length($sequence ));

  print "$header$sequence\n$strand$quality\n";
#  exit;
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
