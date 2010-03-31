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

#my @test_data = ([100,150], [98,101], [87,89], [175,180], [155, 170], [200,250], [185,190]);

#@test_data = ([100,150], [200,250], [300, 350], [180,181], [165, 170],  [285,290], [155,156]);
#@test_data = ([100,150], [50,95], [90,110], [90, 105],  );

my @stack;# = shift @test_data;
#  b_insert(\@stack, \@test_data);

my $dprint = 0;
my $counter = 210;

while(<>) {
  chomp;
  b_insert(\@stack, [[split("\t", $_)]]);
  $dprint++ if ($counter <= 10);
  last if ( ! $counter--);
}


print Dumper(\@stack);
#print "::: ". join( "" ,(map { " (@$_)" } @stack )) ."\n";



# 
# 
# 
# Kim Brugger (30 Mar 2010)
sub b_insert {
  my ( $stack, $input_data ) = @_;

  $$stack[ 0 ][ $SUB ] = [] if ( @$stack == 1 && ! $$stack[ 0 ][ $SUB ] );


  if ( $dprint ){
    print "STACK : ". Dumper( $stack ) ;
    print "INDATA: ".Dumper( $input_data ) ;
#    sleep 1;
  }
  
  foreach my $data ( @$input_data ) {

    $$data[ $SUB ] = [] if ( ! $$data[ $SUB ]);

    # If the stack is empty, just push the first element in and lets get on with it.
    if ( !$stack || !@$stack  ) {
#      print "Up front push \n";
      push @$stack, $data;
      next;
    }
  
    # find the middle point of the array
    
    my ( $left, $right ) = (0, int(@$stack));
    my $middle = floor(($right - $left)/2);
    
#  print "MID :: $middle ( $left, $right)\n";
  print " --> $$data[0], $$data[1]\n" if ( $dprint);
#  print Dumper( \@$stack );
#    print "::: ". join( "" ,(map { " (@$_)" } @$stack )) ."\n";

    $| = 1;
    
    while (1) {
      
      print " $left --> $right || $middle ($$data[0],$$data[1])\n" if ( $dprint);
      print " $$data[ $START ] $$data[ $END ] || $$stack[ $middle ][$START] $$stack[ $middle ][$END]  \n" if ( $dprint);
      sleep 1 if ( $dprint);
      
      if ( $$data[ $END ] < $$stack[ $middle ][$START] ) {
	if ( $middle == 0 ) {
	  unshift @$stack, $data;
	  last;
	}
	elsif ( $$data[ $START ] > $$stack[ $middle - 1 ][$END] ) {
	  splice @$stack, $middle, 0, $data;
	  last;
	}
	else {
	  $right = $middle;
	  last if ( $right <= $left );
	  $middle = $left + floor(($right - $left)/2);
	}
      }
      elsif ( 0 ||  $$data[ $END   ] > $$stack[ $middle ][$START]  &&
	      $$data[ $END   ] < $$stack[ $middle ][$END]  ) {
	print "JOIN TWO BLOCKS 1\n";
#	$dprint = 1;	

	if ( $$data[ $START ] > $$stack[ $middle ][$START] ) {
	  b_insert($stack[ $middle ][$SUB], [$data]);
	  last;
	}
	
	if ( $middle == 0 ||  
	     $$data[ $START ] > $$stack[ $middle - 1][$END] ) {
	  
	  b_insert($$stack[ $middle ][$SUB], [[ $$stack[ $middle ][$START], $$data[ $END ]]]);
	  $$stack[ $middle ][$START] = $$data[ $START ];
	  last;
	}
	elsif ( $middle > 0 && 
		$$data[ $START ] < $$stack[ $middle - 1][$END] ) {
	  b_insert($$stack[ $middle ][$SUB], [[$$data[ $START ], $$stack[ $middle - 1 ][$END]]]);
	  b_insert($$stack[ $middle ][$SUB], [[$$stack[ $middle ][$START], $$data[ $END ]]]);
	  $$stack[ $middle ][$START] = $$stack[ $middle - 1 ][$START];
	  splice @$stack, $middle - 1, 1;
	  last;
	}
	
	
      }
      elsif ( $$data[ $START  ] == $$stack[ $middle ][$START]  &&
	      $$data[ $END    ] == $$stack[ $middle ][$END] ) {

	#Same block
	b_insert( $$stack[ $middle ][$SUB], [[$$data[ $START ], $$stack[ $middle ][$END]]]);
	last;

      }
      elsif ( $$data[ $START   ] > $$stack[ $middle ][$START]  &&
	      $$data[ $START   ] < $$stack[ $middle ][$END]  ) {
	print "JOIN TWO BLOCKS 2\n";

	if ( $$data[ $END   ] <  $$stack[ $middle ][$END] ) {
	  die "Do not go here\n";
	  b_insert($stack[ $middle ][$SUB], [$data]);
	  last;
	}
	
	if ( $middle == @$stack ||  
	     $$data[ $START ] > $$stack[ $middle - 1][$END] ) {
	  
	  b_insert( $$stack[ $middle ][$SUB], [[$$data[ $START ], $$stack[ $middle ][$END]]]);
	  $$stack[ $middle ][$END] = $$data[ $END ];
	  last;
	}
	elsif ( $$data[ $START ] < $$stack[ $middle - 1][$END] ) {
	  b_insert($$stack[ $middle ][$SUB], [[$$data[ $START ], $$stack[ $middle - 1 ][$END]]]);
	  b_insert($$stack[ $middle ][$SUB], [[$$stack[ $middle ][$START], $$data[ $END ]]]);
	  $$stack[ $middle ][$START] = $$stack[ $middle - 1 ][$START];
	  splice @$stack, $middle - 1, 1;
	  last;
	}
	
	
      }
      elsif ($$data[ $START ] > $$stack[ $middle ][$END] ) {
	if ( $middle + 1 == @$stack ) {
	  push @$stack, $data;
	  last;
	}
	elsif ( $$data[ $END ] < $$stack [ $middle + 1 ][ $START ]) {
	  splice @$stack, $middle + 1, 0, $data;
	  last;
	}
	else {
	  $left = $middle;
	  last if ( $right <= $left );
	  $middle = $left + floor(($right - $left)/2);
	}
      }
      
      last if ( $middle > @$stack );
    }
    
  }
}

