#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (16 Nov 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 1;
  if ( $DYNAMIC_LIB_PATHS ) {    
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules" && ! grep /^$path\/modules/, @INC);
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    push @INC, '/home/kb468/easih-toolbox/modules/';
  }
}



use EASIH;
use EASIH::SNPs;
use EASIH::Git;

EASIH::SNPs->New('dbsnp_134');

my $prev = undef;
my $homology_start = undef;
my $homology_snps  = 0;
while(<>) {
  next if (!/^chr\d/ && !/^\d+/);
  chomp;
  my ($pos, $change, $filter, $score, $depth, $genotype, $rest )= split(/\t/);

  next if ($filter ne "PASS");
  next if ( $prev && $prev eq $pos );
#  print "$pos\t$genotype\n";
  if (!$homology_start &&  $genotype =~ /homo/i ) {
    $homology_start = $pos;
    $homology_snps  = 1;
  }
  elsif( $homology_start  &&  $genotype !~ /homo/i ) {
    $homology_start =~ s/chr//;
    $prev =~ s/chr//;
    my $hstart = EASIH::SNPs::CM(split(":", $homology_start));
    my $hend   = EASIH::SNPs::CM(split(":", $prev));
    if ( $hend ne "NA" && $hstart ne "NA" ) {
      my($chr, $start) = split(":", $homology_start);
      ($chr, my  $end) = split(":", $prev);
#      printf("homology block from $homology_start to $prev. Size: %.2f cM (SNPs: $homology_snps)\n", $hend-$hstart) if ( ($hend-$hstart) > 0.5 && $homology_snps>= 3);
      printf("$chr\t$start\t$end\t%.2f cM\t%d bp\t$homology_snps snps\n", $hend-$hstart+1, $end-$start+1) if ( ($hend-$hstart) > 0.5 && $homology_snps>= 3);
#      printf("homology block from $homology_start to $prev. Size: %.2f cM (SNPs: $homology_snps) $hend-$hstart\n", $hend-$hstart);
    }
    $homology_start = undef;
  }
  elsif( $genotype !~ /homo/i ) {
    $homology_start = undef;
  }
  else {
    $homology_snps++;
  }

  $prev = $pos;
}

if( $homology_start ) {
  $homology_start =~ s/chr//;
  $prev =~ s/chr//;
  my $hstart = EASIH::SNPs::CM(split(":", $homology_start));
  my $hend   = EASIH::SNPs::CM(split(":", $prev));
  if ( $hend ne "NA" && $hstart ne "NA") {
    printf("homology block from $homology_start to $prev. Size: %.2f cM\n", $hend-$hstart);
  }
}
