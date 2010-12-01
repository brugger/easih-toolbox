#!/usr/bin/perl 
# 
# transforms a fasta sequence into colourspace...
# 
# 
# Kim Brugger (13 Oct 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $fasta = shift || die "$0 colour_space_string \n";

my %base2colour = (
  'AA' => '0',
  'AC' => '1',
  'AG' => '2',
  'AT' => '3',
  'CA' => '1',
  'CC' => '0',
  'CG' => '3',
  'CT' => '2',
  'GA' => '2',
  'GC' => '3',
  'GG' => '0',
  'GT' => '1',
  'TA' => '3',
  'TC' => '2',
  'TG' => '1',
  'TT' => '0');

my @letters = split( //, $fasta );
my $first_base = $letters[0];
for( my $i = 1; $i < @letters ; $i++ ){
	
  my $encoding = $base2colour{ $first_base.$letters[$i] };
  $first_base = $letters[ $i ];
  $letters[ $i  ] = $encoding;

}
print join("",@letters)."\n";



