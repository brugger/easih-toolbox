#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Mar 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts("r:g:h", \%opts);

usage() if ( $opts{'h'} );


my $RED     = $opts{'r'} || 15;
my $GREEN   = $opts{'g'} || 30;



while(<>) {

  my $name = $_;
  my $seq = <>;
  my $strand  = <>;
  my $qual = <>;

  my @seq = split("", $seq);
  my @qual = split("", $qual);
  
  for(my $i=0;$i<@seq;$i++) {
    my $QV = ord($qual[$i]) - 33;
    if ( $QV < $RED ) {
      $qual[$i] = red($qual[$i]);
      $seq[$i]  = red($seq[$i]);
    }
    elsif ( $QV > $GREEN ) {
      $qual[$i] = green($qual[$i]);
      $seq[$i]  = green($seq[$i]);
    }
    else {
      $qual[$i] = yellow($qual[$i]);
      $seq[$i]  = yellow($seq[$i]);
    }

    $seq  = join("", @seq );
    $qual = join("", @qual);
    
  }
  print "$name$seq$strand$qual";
}



# 
# 
# 
# Kim Brugger (30 Mar 2012)
sub red {
  my ($b) = @_;
  return "\e[0;31m$b\e[0;37m";
}


# 
# 
# 
# Kim Brugger (30 Mar 2012)
sub yellow {
  my ($b) = @_;
  return "\e[0;33m$b\e[0;37m";
}


# 
# 
# 
# Kim Brugger (30 Mar 2012)
sub green {
  my ($b) = @_;
  return "\e[0;32m$b\e[0;37m";
}



# 
# 
# 
# Kim Brugger (30 Mar 2012)
sub usage {
  $0 =~ s/.*\///;
  print STDERR "USAGE $0 : colour bases from a fastq file accoding to base quality\n";
  print STDERR "USAGE $0 : -g<reen cutoff default: 30> -r<ed cutoff, default: 10> -h[elp]\n";
  exit -1;
}
