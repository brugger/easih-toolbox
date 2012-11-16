package EASIH::Seq_mutator;

use strict;
use warnings;
use Data::Dumper;

my $mutation_rate  = 0.001; # 1%
my $deletion_rate  = 0.0005; # 1%
my $insertion_rate = 0.0005; # 1%


# 
# 
# 
# Kim Brugger (21 Feb 2012)
sub set_mutatation_rate {
  my ($new) = @_;
  $mutation_rate = $new if ($new =~ /^\d+\z/ || $new =~ /^\d+.\d+\z/);
  return $mutation_rate;
}

# 
# 
# 
# Kim Brugger (21 Feb 2012)
sub set_deletion_rate {
  my ($new) = @_;
  $deletion_rate = $new if ($new =~ /^\d+\z/ || $new =~ /^\d+.\d+\z/);
  return $deletion_rate;
}

# 
# 
# 
# Kim Brugger (21 Feb 2012)
sub set_insertion_rate {
  my ($new) = @_;
  $mutation_rate = $new if ($new =~ /^\d+\z/ || $new =~ /^\d+.\d+\z/);
  return $insertion_rate;
}




# 
# perl -e '@s=('A','C','G','T');$i = 50; while($i-->0){ print $s[int(rand(4))]} print "\n"'
# 
# Kim Brugger (26 Dec 2011)
sub random_dna {
  my ($bases) = @_;
  $bases ||= 55;
  my @s=('A','C','G','T');
  my $seq = "";
  while($bases-->0){
    $seq .= $s[int(rand(4))]
  }
  
  return $seq;
}


# 
# 
# 
# Kim Brugger (27 Oct 2011)
sub mutate_sequence {
  my ($sequence, $count) = @_;

  my @reads;
  
  $count ||= 2000;
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
    push @reads, $new_read;
  }

  return @reads if ( wantarray );
  return \@reads;
}
