#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (27 Oct 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts("m:d:i:n:h", \%opts);


my $mutation_rate  = 0.01; # 1%
my $deletion_rate  = 0.02; # 1%
my $insertion_rate = 0.06; # 1%


# 
# 
# 
# Kim Brugger (27 Oct 2011)
sub mutate_sequence {
  my ($sequence, $count) = @_;
  
  
  $count ||= 2001;
  my $prefix ||= int(rand(100));
  
  while( $count-- > 0 ) {
    
    my @bases = split(//, $sequence);
    my $new_read;
    
    for( my $i = 0; $i < @bases; $i++) {
      
      my $base = $bases[ $i ] ;
      
      if (rand() < $mutation_rate) {
	while(1) {
	  my $new_base = ('A', 'C', 'G', 'T')[int(rand(3))];
	  
	  if ($new_base  ne $base ) {
	    $base = $new_base;
	    last;
	  }
	}
      }
      elsif (rand() < $deletion_rate) {
	next;
	
      }
      elsif (rand() < $insertion_rate) {
	$base = $base x int(rand(3));
      }
      $new_read .= $base;
    }
    
    print ">$prefix\_$count\_$mutation_rate\_$deletion_rate\_$insertion_rate\n$new_read\n";
    
  }
}
