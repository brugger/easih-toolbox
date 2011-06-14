#!/usr/bin/perl 
# 
# Take a colour space line and transform it into base space.
# 
# 
# Kim Brugger (10 Jun 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


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

#die "program csfasta shift\nif shift is 1, the first base is omitted in the output\n" if (@ARGV < 1);

my $csfasta = shift || die "$0 colour_space_string or csfasta file\n";

if (! -e $csfasta ) {
  print decode( $csfasta) ."\n";
}
else {
  my ($name, $seq);
  open(my $in, $csfasta) || die "Could not open 'csfasta': $!\n";
  while(<$in>) {
    
#    print ;

    if (/^\>/) {
      $name = $_;
      next;
    }
    elsif (/^[ACGT]\d+\Z/ && $name) {
      print $name . decode($_) . "\n";
      $name = undef;
    }
  }
}    


  
# 
# 
# 
# Kim Brugger (29 Mar 2011)
sub decode {
  my ($cs) = @_;

  $cs =~ s/\n//g;

  my @letters = split( //, $cs );
  my $first_base = $letters[0];
  for( my $i = 1; $i < @letters ; $i++ ){
    
    my $colour = $letters[$i];
    my $encoding = $first_base.$colour;
    $first_base = $colourspace{ $encoding };
    $letters[ $i ] = $first_base;    
  }
  shift @letters;
  return join("",@letters);
}
  
  


