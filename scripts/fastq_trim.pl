#!/usr/bin/perl 
# 
# Trims fq entries to the requested length.
# 
# 
# Kim Brugger (16 May 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::Git;

my %opts;
getopts("l:hi:o:", \%opts);

my $infile   = $opts{i} || "";;
my $outfile  = $opts{o} || "";;
my $length   = $opts{l} || Usage();

open (*STDIN,  $infile) || die "Could not open '$infile': $!\n" if ( $infile );
open (*STDOUT, ">$outfile") || die "Could not write to '$outfile': $!\n" if ( $outfile );

while(my $header = <> ) {
  my $sequence   = <>;
  my $strand     = <>;
  my $quality    = <>;

  die "fastq file is either not a fq-file or contains extra lines.\n"  
      if ( $header !~ /^\@/ || $strand !~ /^[+-]/ );
  

  die "Sequence length is shorter than requested length\n" if ( length($sequence) < $length);

  $sequence = substr( $sequence, 0, $length);
  $quality  = substr( $quality,  0, $length);


  print "$header$sequence\n$strand$quality\n";
#  exit;
}



# 
# 
# 
# Kim Brugger (16 May 2011)
sub Usage {

  $0 =~ s/.*\///;
  print "USAGE: $0 trims the entries of a fastq to a requested length\n";
  print "USAGE: $0 -l[ength wanted] -i<nfile, default stdin> -o<utfile, default stdout>\n";
  exit -1;
  
}
