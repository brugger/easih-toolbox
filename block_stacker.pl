#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Mar 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
no warnings 'recursion';
use Data::Dumper;
use POSIX qw(ceil floor);


my $START = 0;
my $END   = 1;
my $SUB   = 2;


my $VERBOSE_LEVEL = 00;

my @test_data = ([100,150], [98,101], [87,89], [175,180], [155, 170], [200,250], [185,190]);

@test_data = ([100,150], [200,250], [300, 350], [180,181], [165, 170],  [285,290], [155,156]);
@test_data = ([100,150], [50,95], [90,110], [90, 105],  );

my @stack;# = shift @test_data;

#my $stack;

my $dprint = 2;
#my $counter = 520;

use Time::HiRes;

my $t0 = Time::HiRes::gettimeofday();

print "-"x60 . "\n";
print " "x25 . "START RUN\n";
print "-"x60 . "\n";


if (1) {
  my $counter = 0;
  while(<>) {
    chomp;
    print "STATUS: $counter entries, array length: " . @stack ." \n" if ( $counter++ % 10000 == 0);
    
    b_insert(\@stack, [split("\t", $_)]);
#    print_stack( \@stack, 0 );
#    print "--------------------------\n\n";
 #   $dprint++ if ($counter >= 60);
#    $VERBOSE_LEVEL = 10 if ($counter >= 1051);
#    last if ( ! $counter--);
  }
}
else {
  b_inserts(\@stack, \@test_data);
}

my $t1 = Time::HiRes::gettimeofday(); 

print "-"x60 . "\n";
print " "x25 . "ENDED RUN\n";
print " "x20 . "RUNTIME :: ". int(( $t1 - $t0)/60 ) . " min\n";
print "-"x60 . "\n";
print_stack( \@stack, 100 );

my($nodes, $max_depth, $troughs) = stack_stats( \@stack );

print "Stack stats :: nodes= $nodes, max depth= $max_depth, troughs= $troughs \n";

# 
# 
# 
# Kim Brugger (01 Apr 2010)
sub stack_stats {
  my ( $stack_ref, $level, $max_depth, $troughs  ) = @_;

  my $nodes = 0;
  $max_depth ||= 1;
  $level     ||= 1;
  $troughs   ||= 0;

  $max_depth = $level if ( $level > $max_depth);

  foreach my $field ( @$stack_ref ) {    
    my $cutoff = 100;
    next if ( $$field[ $END ] - $$field[ $START ] + 1 < $cutoff);
    
    $troughs++ if ( $level == 1);

    $nodes++;

    if ( @{$$field[ $SUB ]} > 0 ) {
      my ($sub_nodes, $l_max_depth) = stack_stats($$field[ $SUB ], $level + 1, $max_depth);

      $nodes += $sub_nodes;
      $max_depth = $l_max_depth if ( $max_depth < $l_max_depth);
    }
    
  }
  
  $troughs--;

  return( $nodes, $max_depth, $troughs);
}


# 
# 
# 
# Kim Brugger (31 Mar 2010)
sub print_stack {
  my ( $stack_ref, $cutoff, $level ) = @_;

  $level ||= 1;

  foreach my $field ( @$stack_ref ) {    

    next if ( $$field[ $END ] - $$field[ $START ] + 1 < $cutoff);

    print "|"x$level . "-$$field[ $START ]->$$field[ $END ] (".($$field[ $END ] - $$field[ $START ] + 1).")\n";
    if ( @{$$field[ $SUB ]} > 0 ) {
      print "|"x$level."\\\n";
      print_stack( $$field[ $SUB ], $cutoff, $level + 1);
      print "|"x$level."/\n";
    }
  }
  
}



# 
# Simple wrapper for inserting multiple data fields.
# 
# Kim Brugger (31 Mar 2010)
sub b_inserts {
  my ( $stack, $input_data ) = @_;
  
  foreach my $data ( @$input_data ) {
    $stack = b_insert( $stack, $data);
  }
  
}



# 
# 
# 
# Kim Brugger (30 Mar 2010)
sub b_insert {
  my ( $stack_ref, $data ) = @_;

  if ( 1 && $dprint ){
#    verbose("STACK_REF : ". Dumper( $stack_ref ), 2);
#    verbose("INDATA: ". Dumper( $data ), 2);
  }

 RERUN:

#  $$stack_ref[ 0 ][ $SUB ] = [] if ( @$stack_ref == 1 && ! $$stack_ref[ 0 ][ $SUB ] );


  

  # initialise the sub array value so we can use it as a reference later, otherwise the 
  # information is lost.
  $$data[ $SUB ] = [] if ( ! $$data[ $SUB ]);


  # If the stack_ref is empty, just push the first element in and lets get on with it.
  if ( !$stack_ref || !@$stack_ref  ) {
    verbose("Up front push \n", 2);
    push @$stack_ref, $data;
#    print Dumper( $stack_ref );
#    print Dumper( \@stack );
    return $stack_ref;
  }
  
  # set the start and end of the array and find the 
  # the middle of the array
  my ( $left, $right ) = (0, int(@$stack_ref));
  my $middle = floor(($right - $left)/2);
    
#  verbose("MIDDLE $middle ( $left, $right)\n", 1) if ( $dprint);
#  verbose("BLOCK  $$data[0], $$data[0]\n", 1) if ( $dprint);
#    print "::: ". join( "" ,(map { " (@$_)" } @$stack_ref )) ."\n";

  # Flush the buffer constantly
  $| = 1;

  my $loop_counter = 0;
    
  while (1) {

    if ( $loop_counter++ > @$stack_ref ) {
      $VERBOSE_LEVEL = 10;
      verbose("MIDDLE $middle ( $left, $right)\n", 1);
      verbose(" $$data[ $START ] $$data[ $END ] <=> $$stack_ref[ $middle ][$START] $$stack_ref[ $middle ][$END]\n", 1);
#    sleep 1 if ( $dprint);
      exit;
    }

    verbose("MIDDLE $middle ( $left, $right)\n", 1);
    verbose(" $$data[ $START ] $$data[ $END ] <=> $$stack_ref[ $middle ][$START] $$stack_ref[ $middle ][$END]\n", 1);

    
    # The new block is to the left of the middle.
    if ( $$data[ $END ] < $$stack_ref[ $middle ][$START] ) {
      # there a no blocks further upstream, so just unshift it in to the beginning of the array
      if ( $middle == 0 ) {
	unshift @$stack_ref, $data;
	last;
      }
      # the block is between the middle block and the block just downstream
      elsif ( $$data[ $START ] > $$stack_ref[ $middle - 1 ][$END] ) {
	splice @$stack_ref, $middle, 0, $data;
	last;
      }
      # no cigar, make the next jump
      else {
	$right = $middle;
	last if ( $right <= $left );
	$middle = $left + floor(($right - $left)/2);
      }
    }
    # The new block is to the right of the middle.
    elsif ($$data[ $START ] > $$stack_ref[ $middle ][$END] ) {
      # there are no blocks further downstrem, so just push onto the end of the array
      if ( $middle + 1 == @$stack_ref ) {
	push @$stack_ref, $data;
	last;
      }
      # the block is between the middle on the the next one downstream so we will just slot it in.
      elsif ( $$data[ $END ] < $$stack_ref [ $middle + 1 ][ $START ]) {
	splice @$stack_ref, $middle + 1, 0, $data;
	last;
      }
      # No cigar, so make the jump.
      else {
	$left = $middle;
	last if ( $right <= $left );
	$middle = $left + floor(($right - $left)/2);
      }
    }
    #
    # Now things gets interesting, we here start to calculate
    # overlapping and contained regions.
    #
    
    # this is a contained or identical block
    elsif ( $$data[ $START  ] >= $$stack_ref[ $middle ][ $START ]  &&
	    $$data[ $END    ] <= $$stack_ref[ $middle ][ $END   ] ) {
      
      verbose("CONTAINED BLOCK", 2);
      b_insert( $$stack_ref[ $middle ][$SUB], $data);
      last;
    }
    # The new block is a container for the middle block
    elsif ( $$data[ $START  ] <= $$stack_ref[ $middle ][ $START ]  &&
	    $$data[ $END    ] >= $$stack_ref[ $middle ][ $END   ] ) {

#      $VERBOSE_LEVEL = 10;
      verbose("CONTAINER BLOCK", 2);



      # If there is only one element in the array we can loose the reference to our data
      # so the values are just replaced, and things work fine.
      if ( @$stack_ref == 1 ) {

	verbose("SINGLE ENTRY ARRAY \n", 2);
#	print "REFE: " . Dumper( $stack_ref );
#	print "DATA: " . Dumper( $data );

	my (@middle_data) = ($$stack_ref[0][$START], $$stack_ref[0][$END], $$stack_ref[0][$SUB]);

#	($$stack_ref[ $START ], $$stack_ref[ $END ], $$stack_ref[ $SUB ]) = ($$data[0], $$data[1], $$data[2]);
#	@$stack_ref = @$data;
	
	($$stack_ref[0][$START], $$stack_ref[0][$END], $$stack_ref[0][$SUB]) = ($$data[$START], $$data[$END], $$data[$SUB]);
	
#	print "DATA: " .Dumper( $data );
#	print "REFE: " . Dumper( $stack_ref );
#	print "ORIN: " . Dumper( \@stack );
#	print "MIDD: " .Dumper( \@middle_data );
#	$dprint++;
	b_insert( $stack_ref, [@middle_data]);
#	exit;
	return;
      }

#      $VERBOSE_LEVEL = 10;
      verbose("MULTIPLE ENTRY ARRAY \n", 2);
      # pick the middle block out of the array and insert it into the current block;
      my $middle_block = splice @$stack_ref, $middle, 1;
#      print "CONTAINED : " . Dumper( $middle_block);
#      print "PRE-DATA : " . Dumper( $data);
#      print "STACK    : " . Dumper( $stack_ref );
      b_insert([ $data], $middle_block);
#      print "POST-DATA : " . Dumper( $data);


      #As this block could overlap with other blocks, rerun the loop
      verbose("RERUNNING LOOP\n", 2);
      goto RERUN;
     
      
      next;
    }
    # New block overlaps to the left of the middle block
    elsif ( $$data[ $END   ] >= $$stack_ref[ $middle ][ $START ]  &&
	    $$data[ $END   ] <= $$stack_ref[ $middle ][ $END   ]  ) {
      verbose("JOIN TWO BLOCKS LEFT\n", 2);
      
      # We are either at the start of the array or 
      # the block does not overlap with the next block downstream
      if ( $middle == 0 ||  
	   $$data[ $START ] >= $$stack_ref[ $middle - 1][ $END ] ) {

	verbose("SINGLE JOIN\n", 1);
	b_insert($$stack_ref[ $middle ][$SUB], [ $$stack_ref[ $middle ][$START], $$data[ $END ]]);
	$$stack_ref[ $middle ][ $START ] = $$data[ $START ];
	last;
      }
      # Not at the start of the array, and we are overlapping with the next block downstream
      # but the downstream block is _not_ contained in this new block
      elsif ( $middle > 0 && 
	      $$data[ $START ] <= $$stack_ref[ $middle - 1][ $END   ] && 
	      $$data[ $START ] >= $$stack_ref[ $middle - 1][ $START ] ) {

	verbose("DOUBLE JOIN\n", 1);
	# insert the two new overlaps
	b_insert($$stack_ref[ $middle ][$SUB], [$$data[ $START ], $$stack_ref[ $middle - 1 ][$END]]);
	b_insert($$stack_ref[ $middle ][$SUB], [$$stack_ref[ $middle ][$START], $$data[ $END ]]);
	# adjust the start position of the block
	$$stack_ref[ $middle ][ $START ] = $$stack_ref[ $middle - 1 ][ $START ];
	# and remove the block downstream
	splice @$stack_ref, $middle - 1, 1;
	last;
      }
    }
    # New block overlaps to the right of the middle block
    elsif ( $$data[ $START   ] >= $$stack_ref[ $middle ][$START]  &&
	    $$data[ $START   ] <= $$stack_ref[ $middle ][$END]  ) {
      verbose("JOIN TWO BLOCKS RIGHT\n", 2);
      
      if ( $$data[ $END   ] <  $$stack_ref[ $middle ][$END] ) {
	die "Do not go here\n";
	b_insert($$stack_ref[ $middle ][$SUB], [$data]);
	last;
      }

      # 
      if ( $middle + 1 == @$stack_ref ||  
	   $$data[ $START ] >= $$stack_ref[ $middle - 1][$END] ) {
	
	verbose("SINGLE JOIN\n", 1);
	b_insert( $$stack_ref[ $middle ][$SUB], [$$data[ $START ], $$stack_ref[ $middle ][$END]]);
	$$stack_ref[ $middle ][$END] = $$data[ $END ];
	last;
      }
      elsif ( $$data[ $START ] < $$stack_ref[ $middle - 1][$END] ) {
	verbose("DOUBLE JOIN\n", 1);
	b_insert($$stack_ref[ $middle ][$SUB], [$$data[ $START ], $$stack_ref[ $middle - 1 ][$END]]);
	b_insert($$stack_ref[ $middle ][$SUB], [$$stack_ref[ $middle ][$START], $$data[ $END ]]);
	$$stack_ref[ $middle ][$START] = $$stack_ref[ $middle - 1 ][$START];
	splice @$stack_ref, $middle - 1, 1;
	last;
      }
      
      
    }
    
    last if ( $middle > @$stack_ref );
  }
    
  return;
}




# 
# 
# 
# Kim Brugger (31 Mar 2010)
sub verbose {
  my ($message, $level) = @_;
  return if ( $level > $VERBOSE_LEVEL);
  $message =~ s/\n+\Z//g;
  print "MESS " . ":"x$level . " $message\n";
#  sleep 1;
}

