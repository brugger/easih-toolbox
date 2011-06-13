#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (17 Mar 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts("m:c:h", \%opts);


usage() if ( $opts{h});

my $max_mismatches = $opts{'m'} || 10;
my $cluster_min    = $opts{'c'} || 100;

my %s = ( A => '0', 
	  0 => 'A',
	  C => '1', 
	  1 => 'C',
	  G => '2', 
	  2 => 'G',
	  T => '3',
	  3 => 'T',
	  N => '4');

my %reads;
# for debugging
my $limit = 40000;

while( <> ) {
  my $name       = $_;
  my $sequence   = <>;
  my $strand     = <>;
  my $quality    = <>;
  
  chomp($name);
  chomp($sequence );
  chomp($quality);

  push @{$reads{ $sequence }}, [$name, $quality];
#  last if ( ! $limit--);
}


# Cluster the similar reads together.
my (%merged);
foreach my $seq ( sort { @{$reads{$b}} <=> @{$reads{$a}} } keys %reads) {

  next if ( $seq =~ /^\.+\Z/);
  next if ( $seq =~ /\./);
  
  if ( ! keys %merged ) {
    push @{$merged{ $seq }}, $seq;
  }
  else {
    my $seq2 = similar_seqs( $seq, keys %merged);
    if ( $seq2 ) {
      push @{$merged{ $seq2 }}, $seq;
      next;
    }
    
    push @{$merged{ $seq }}, $seq;
  }
}



foreach my $main_seq ( keys %merged ) {

  my $clustered_reads = 0;
  map { $clustered_reads += @{$reads{$_}}} @{$merged{$main_seq}};

  next if ($clustered_reads < $cluster_min );



  my @data;
  # make a variation profile for the cluster
  foreach my $seq ( @{$merged{$main_seq}} ) {
    my $count = @{$reads{ $seq }};
    my @seq   = split("", $seq);

    for(my $i=0;$i<@seq;$i++) {
      $data[$i][$s{$seq[$i]}] += $count;
    }
  }

  # Calculate the regions where corrections can be made.
  my @corrections;
  for(my $i=0;$i<@data;$i++) {
    
    my $d =  $data[$i];
    my ($a1, $a2) = major( $d );
    
    if ( $a1 && $a2 ) {
      $corrections[$i] = [$a1, $a2];
    }
    elsif ( $a1 ) {
      $corrections[$i] = [$a1];
    }
    
  }

  # and do the corrections.
  foreach my $seq ( @{$merged{$main_seq}} ) {
    my @seq = split("", $seq);
    
    foreach my $read ( @{$reads{$seq}} ) {
      my $corrected_seq = "";
      my @qual = split("", $$read[1]);
      for(my $i=0;$i<@seq;$i++) {
	if ( @{$corrections[$i]} == 1 ) {
	  if ( $corrections[$i][0] eq $seq[$i]) {
	    $corrected_seq .= $seq[$i];
	  }
	  else {
	    $corrected_seq .= lc($corrections[$i][0]);
	    $qual[$i] = "!"
	  }
	}
	else {
	  $corrected_seq .= $seq[$i];
	}
      }
      
#      print "$$read[0]\n$corrected_seq\n+\n".join("", @qual)."\n";
#      print "$$read[0]\n$corrected_seq\n";
      print "$corrected_seq\n";
    }
  }
#  exit;
}




# 
# 
# 
# Kim Brugger (17 Mar 2011)
sub major {
  my ( $d ) = @_;

  $$d[ $s{A}] ||= 0;
  $$d[ $s{C}] ||= 0;
  $$d[ $s{G}] ||= 0;
  $$d[ $s{T}] ||= 0;

  my $A = $$d[ $s{A}];
  my $C = $$d[ $s{C}];
  my $G = $$d[ $s{G}];
  my $T = $$d[ $s{T}];

  my $max = $A + $C + $G + $T;
  $max = $max - $max*0.02;

  if ( $A > $max ) {
#    print "A: $A\n";
    return("A");
  }    
  elsif ( $C > $max ) {
#    print "C: $C\n";
    return("C");
  }    
  elsif ( $G > $max ) {
#    print "G: $G\n";
    return("G");
  }    
  elsif ( $T > $max ) {
#    print "T: $T\n";
    return("T");
  }    
  else {

    my @r = sort{ $$d[$b] <=> $$d[$a] } (0..3);
    my ($first, $second) = ($r[0], $r[1]);
#    print "$s{$first}$s{$second}: $$d[$first]/$$d[$second]\n";
    return ($s{$first},$s{$second});
    
  }
  
  return (undef, undef);
  
}

  

# 
# 
# 
# Kim Brugger (06 Jan 2011)
sub similar_seqs {
  my ($seq1, @seq2s) = @_;

  foreach my $seq2 ( @seq2s ) {

    my $mask = $seq1 ^ $seq2;
    my $diffs = 0;
    while ($mask =~ /[^\0]/g) {
      $diffs++;
    }    

    return $seq2 if ( $diffs <= $max_mismatches);
  }

  return undef;
}


# 
# 
# 
# Kim Brugger (18 Mar 2011)
sub usage {
  $0 =~ s/.*\///;
  print "SAET for fastq files, similar to what can be done with SOLiD data\n";
  print "USAGE: $0 -m[ax ismatches for two sequences to be clustered together (default: 10)] -c[luster min size (default: 100)]\n";
  print "Devised and implemented by the now soiled: Kim Brugger (18 Mar 2011), contact: kim.brugger\@easih.ac.uk\n";
  exit -1;
}
