#!/usr/bin/perl 
# 
# Take a colour space line and transform it into base space.
# 
# 
# Kim Brugger (10 Jun 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;



my $csfasta = shift || die "$0 colour_space_string \n";
my %colourspace = (
    "A0" => "A",
    "C0" => "C",
    "G0" => "G",
    "T0" => "T",
    "A1" => "C",
    "C1" => "A",
    "G2" => "A",
    "A2" => "G",
    "A3" => "T",
    "T3" => "A",
    "C2" => "T",
    "T2" => "C",
    "C3" => "G",
    "G3" => "C",
    "G1" => "T",
    "T1" => "G" );

my @letters = split( //, $csfasta );
my $first_base = $letters[0];
for( my $i = 1; $i < @letters ; $i++ ){
	
  my $colour = $letters[$i];
  my $encoding = $first_base.$colour;
  $first_base = $colourspace{ $encoding };
  $letters[ $i ] = $first_base;    
}
shift @letters;
print join("",@letters)."\n";



