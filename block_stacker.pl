#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Mar 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use POSIX qw(ceil floor);


my $START = 0;
my $END   = 0;
my $LEVEL = 0;

my @data = ([100,150], [98,101], [87,89], [175,180], [155, 170], [200,250], [185,190]);

@data = ([100,150], [200,250], [300, 350], [180,181], [165, 170],  [285,290], [155,156]);

my @stack = shift @data;

foreach my $data ( @data ) {


  # this is the first element, so just push it in.
  if ( !@stack ) {
    push @stack, $data;
    next;
  }
  
  # find the middle point of the array

  my ( $left, $right ) = (0, int(@stack));
  my $middle = floor(($right - $left)/2);

#  print "MID :: $middle ( $left, $right)\n";
#  print " --> $$data[0], $$data[1]\n";
#  print Dumper( \@stack );
  print "::: ". join( "" ,(map { " (@$_)" } @stack )) ."\n";
  
  while (1) {

    print " $left --> $right || $middle ($$data[0],$$data[1])\n";
    
    if ( $$data[ $END ] < $stack[ $middle ][$START] ) {
      if ( $middle == 0 ) {
	unshift @stack, $data;
	last;
      }
      elsif ( $$data[ $START ] > $stack[ $middle - 1 ][$END] ) {
	splice @stack, $middle, 0, $data;
	last;
      }
      else {
	$right = $middle;
	last if ( $right <= $left );
	$middle = $left + floor(($right - $left)/2);
      }
    }
    elsif ($$data[ $START ] > $stack[ $middle ][$END] ) {
      if ( $middle + 1 == @stack ) {
	push @stack, $data;
	last;
      }
      elsif ( $$data[ $END ] < $stack [ $middle + 1 ][ $START ]) {
	splice @stack, $middle + 1, 0, $data;
	last;
	      }
      else {
	$left = $middle;
	last if ( $right <= $left );
	$middle = $left + floor(($right - $left)/2);
      }
    }
  
    last if ( $middle > @stack );
  }
  
  
}


#print Dumper(\@stack);
print "::: ". join( "" ,(map { " (@$_)" } @stack )) ."\n";
