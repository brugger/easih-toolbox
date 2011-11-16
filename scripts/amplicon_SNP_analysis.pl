#!/usr/bin/perl -w
# 
# 
# 
# 
# Kim Brugger (02 Aug 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::ChiSquare;
use EASIH::Fastq;
my $infile = shift;

my %res;


my $in = EASIH::Fastq::open( $infile );

while( 1 ) {
  my (undef, $seq, undef, undef) = EASIH::Fastq::next($in);
  last if ( ! $seq );
  chomp( $seq );
  $res{ $seq }++;
}

# fetch the two largest sequence populations
my @best = (sort { $res{ $b } <=> $res{ $a }} keys %res)[0..1];

my ($count1, $count2) = ($res{ $best[0]}, $res{ $best[1]});

my $lcount1 = 10*log($count1);
my $lcount2 = 10*log($count2);

my  ($chisquare, $degrees_of_freedom, $chip) =  EASIH::ChiSquare::chisquare_raw( [$lcount1, $lcount2] );


print "SNP is Heterozygous\tP=$chip,read1=$count1,read2=$count2\n" if ( $chip*100 >= 5);
print "SNP is Homozygous\tP=$chip,read1=$count1,read2=$count2\n" if ( $chip*100 < 5);


