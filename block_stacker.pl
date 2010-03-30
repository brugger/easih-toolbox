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
my $END   = 1;
my $SUB   = 2;

my @data = ([100,150], [98,101], [87,89], [175,180], [155, 170], [200,250], [185,190]);

@data = ([100,150], [200,250], [300, 350], [180,181], [165, 170],  [285,290], [155,156]);
@data = ([100,150], [50,95], [90,110], [91, 105] );

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
    print "    if ( $$data[ $END ] < $stack[ $middle ][$START] ) { \n";

    
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
    elsif ( $$data[ $END   ] > $stack[ $middle ][$START] ) {
      print "JOIN TWO BLOCKS \n";
      
      if ( $$data[ $START ]>  $stack[ $middle ][$START] && 
	   $$data[ $END   ] > $stack[ $middle ][$START] ) {
	push @{$stack[ $middle ][$SUB]}, $data;
	last;
      }

      if ( $middle == 0 ||  
	   $$data[ $START ] > $stack[ $middle - 1][$END] ) {
	
	push @{$stack[ $middle ][$SUB]}, [ $stack[ $middle ][$START], $$data[ $END ]];
	$stack[ $middle ][$START] = $$data[ $START ];
	last;
      }
      elsif ( $middle > 0 && 
	   $$data[ $START ] < $stack[ $middle - 1][$END] ) {
	push @{$stack[ $middle ][$SUB]}, [$$data[ $START ], $stack[ $middle - 1 ][$END]];
	push @{$stack[ $middle ][$SUB]}, [$stack[ $middle ][$START], $$data[ $END ]];
	$stack[ $middle ][$START] = $stack[ $middle - 1 ][$START];
	print 
	splice @stack, $middle - 1, 1;
	last;
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


print Dumper(\@stack);
#print "::: ". join( "" ,(map { " (@$_)" } @stack )) ."\n";
