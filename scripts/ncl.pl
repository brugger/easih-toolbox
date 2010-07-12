#!/usr/bin/perl


# Nested Containment List (NCList): A new algorithm for accelerating
# interval query of genome alignment and interval databases
# http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btl647v1


# Nested lists are composed of intervals defined by begin and end positions,
# and any interval that is contained within another interval is rooted to this.
# Thus, fast interval lookups can be performed using binary search.


use warnings;
use strict;
use Data::Dumper;
use Storable;
no warnings 'recursion';


my ( $features, $nc_list );


use Time::HiRes;

my $t0 = Time::HiRes::gettimeofday();


my $counter = 0;
while(<>) {
  chomp;
  
  next if (/^\z/);
  push @$features, [split("\t", $_)];
}



$nc_list = nc_list_create( $features, 0, 1, 2 );

my $t1 = Time::HiRes::gettimeofday(); 
print "RUNTIME :: ". int(( $t1 - $t0) ) . " sek(s)\n";

print Dumper( $nc_list );


exit;





# Martin A. Hansen, February 2010.
# Creates a Nested Containment (NC) list from a stack of features.
# The features consits of an AoA where beg and end specifies the
# elements containing the begin and end position of a feature, and
# index specifies a element used for nesting lists.
sub nc_list_create {
  my ( $features, # list of features AoA
       $beg, # feature element with begin position
       $end, # feature element with end position
       $index, # feature element with index position
      ) = @_;


  my ( $nc_list );


  @{ $features } = sort { $a->[ $beg ] <=> $b->[ $beg ] or $b->[ $end ] <=> $a->[ $end ] } @{ $features };
  $nc_list = [ shift @{ $features } ];
  map { nc_list_add_ns( $nc_list, $_, $end, $index ) } @{ $features };

  return wantarray ? @{ $nc_list } : $nc_list;
}



# Martin A. Hansen, February 2010.
# Recursively construct a Nested Containment (NC) list by added
# a given feature to an existing NC list.
sub nc_list_add_ns {
  my ( $nc_list, # NC list
       $feat, # feature (AoA)
       $end, # feature element with end position
       $index # feature element with index position
      ) = @_;


  # feature is nested.
  if ( $feat->[ $end ] <= $nc_list->[ -1 ]->[ $end ] ) {
    # sublist exists so recurse to this.
    if ( defined $nc_list->[ -1 ]->[ $index ] ) { 
      nc_list_add( $nc_list->[ -1 ]->[ $index ], $feat, $end, $index );
    } 
    else {
      # creating a new sublist.
      $nc_list->[ -1 ]->[ $index ] = [ $feat ]; 
    }
  }
  else {
    push @{ $nc_list }, $feat;
  }
}



# Martin A. Hansen, February 2010.
# Recursively construct a Nested Containment (NC) list by added
# a given feature to an existing NC list.
sub nc_list_add {
  my ( $nc_list, # NC list
       $feat, # feature (AoA)
       $end, # feature element with end position
       $index # feature element with index position
      ) = @_;

  # feature is nested.
  if ( $feat->[ $end ] <= $nc_list->[ -1 ]->[ $end ] ) {
    # sublist exists so recurse to this.
    if ( defined $nc_list->[ -1 ]->[ $index ] ) { 
      nc_list_add( $nc_list->[ -1 ]->[ $index ], $feat, $end, $index );
    } 
    else {
      # creating a new sublist.
      $nc_list->[ -1 ]->[ $index ] = [ $feat ]; 
    }
  }
  else {
    push @{ $nc_list }, $feat;
  }
}
