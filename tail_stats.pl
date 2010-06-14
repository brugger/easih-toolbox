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

use lib '/home/kb468/easih/modules/';
use EASIH::Profile;

my $profile = EASIH::Profile->New;

my $bam_file = shift || die "No bamfile supplied\n";

$profile->from_bam($bam_file, "chrX");
$profile->import_interest_regions('/home/kb468/100_genes_plus_conserved_regions.txt');
my ($off_target_profile, $on_target_profiles) = $profile->split_by_interest_regions();

my @five_tail;
my @three_tail;

foreach my $otp ( @$on_target_profiles ) {

  if ($$otp{five}) {

    my $profile =  $$otp{five}->to_text();
    for( my $i= 0; $i < @$profile; $i++ ) {
      push @{$five_tail[ $i ]}, $$profile[ $i ][ 1 ];
    }
  }

  if ($$otp{three}) {

    my $profile =  $$otp{three}->to_text();
    for( my $i= 0; $i < @$profile; $i++ ) {
      push @{$three_tail[ $i ]}, $$profile[ $i ][ 1 ];
    }
  }
}

for( my $i = 0; $i < @five_tail; $i++ ) {
  print "-".(@five_tail - $i)."\t". int((eval join '+', @{$five_tail[ $i ]})/@{$five_tail[ $i ]}) . "\n";;
}

for( my $i = 0; $i < @three_tail; $i++ ) {
  print "".($i+1)."\t". int((eval join '+', @{$three_tail[ $i ]})/@{$three_tail[ $i ]}) . "\n";;
}
