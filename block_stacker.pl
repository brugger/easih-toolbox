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

my @data = ([100,150], [90,91], [95,97], [92,94], [87,89], [98,99], [175,180], [155, 170], [200,250]);

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

  print "MID :: $middle ( $left, $right)\n";
  print " --> $$data[0], $$data[1]\n";
#  print Dumper( \@stack );
  print "ARRAY:: ". join( "" ,(map { " (@$_)" } @stack )) ."\n";
  
  while (1) {
    
    if ( $$data[ $END ] < $stack[ $middle ][$START] ) {
      if ( $middle == 0 ) {
#	print "UNSHIFT\n";
	unshift @stack, $data;
	last;
      }
      elsif ( $$data[ $START ] > $stack[ $middle - 1 ][$END] ) {
#	print "SPLICE \n";
	splice @stack, $middle, 0, $data;
	last;
      }
      else {
	$right = $middle;
	$middle = $left + floor(($right - $left)/2);
#	print "MID2 :: $middle ($left, $right)\n";

      }
    }
    elsif ($$data[ $START ] > $stack[ $middle ][$END] ) {
      if ( $middle + 1 == @stack ) {
#	print "PUSHING \n";
	push @stack, $data;
	last;
      }
      else {
	$left = $middle;
	$middle = $left + floor(($right - $left)/2);
#	print "MID3 :: $middle ($left, $right)\n";
      }
    }
  
    last if ( $middle > @stack );
  }
  
  
}


#print Dumper(\@stack);
print "SORTED ARRAY:: ". join( "" ,(map { " (@$_)" } @stack )) ."\n";
