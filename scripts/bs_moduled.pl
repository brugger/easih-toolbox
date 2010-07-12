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
use EASIH::Timer;
use EASIH::Stacker;
use EASIH::Profile;

my $profile = EASIH::Profile->New;

if ( 1 ) {
  $profile->from_bam("/home/kb468/solid/0354/XLMR.paired.bioscope.bam", "chrX");
  $profile->import_interest_regions('/home/kb468/100_genes_plus_conserved_regions.txt');
}
else {
  $profile->from_bam("/home/kb468/solid/0354/XLMR.single.bwa.bam", "X");
  $profile->import_interest_regions_gff('/home/kb468/converted_100_genes_plus_conserved_regions.gff');
}
my ($off_target_profile, $on_target_profiles) = $profile->split_by_interest_regions();

my $summed_5_length  = 0;
my $summed_3_length  = 0;
my $summed_on_length = 0;
my $summed_5_mean    = 0;
my $summed_3_mean    = 0;
my $summed_on_mean   = 0;
my $mappings_5       = 0;
my $mappings_3       = 0;
my $mappings_on      = 0;

my $zero_mappings    = 0;

#print Dumper($on_target_profiles );
foreach my $otp ( @$on_target_profiles ) {

#  print Dumper( $otp);
  
  if ($$otp{five}) {
    $mappings_5++;
    $summed_5_length += $$otp{five}->length();
    $summed_5_mean   += $$otp{five}->mean_coverage();
  }

  if ($$otp{on}) {
    my $on_mean = $$otp{on}->mean_coverage();
    if ( $on_mean > 0 ) {
      $mappings_on++;
      $summed_on_length += $$otp{on}->length();
      $summed_on_mean   += $on_mean;
    }
    else {
      $zero_mappings++;
    }
  }

  if ($$otp{three}) {
    $mappings_3++;
    $summed_3_length += $$otp{three}->length();
    $summed_3_mean   += $$otp{three}->mean_coverage();
  }

  if (0) {
    print "region: " . join('->', $$otp{on}->range) . "\n";
    print "5' mean: " . $$otp{five}->mean_coverage() . "\n"  if ($$otp{five});
    print "5' length: " . $$otp{five}->length() . "\n"  if ($$otp{five});
    print "on mean: " . $$otp{on}->mean_coverage() . "\n";
    print "3' mean: " . $$otp{three}->mean_coverage() . "\n" if ($$otp{three});
    print "3' length: " . $$otp{three}->length() . "\n"  if ($$otp{five});
  }
}

print "mapping stats: \n";
printf("on: mappings    : %7d\n",  $mappings_on) if ($mappings_on);
printf("on:avg length   : %7.2f\n", $summed_on_length/ $mappings_on) if ($mappings_on);
printf("on:avg coverage : %7.2f\n", $summed_on_mean/ $mappings_on)   if ($mappings_on);
printf("5 : mappings    : %7d\n",  $mappings_5) if ($mappings_on);
printf("5 :avg length   : %7.2f\n", $summed_5_length/ $mappings_5)   if ($mappings_5);
printf("5 :avg coverage : %7.2f\n", $summed_5_mean/ $mappings_5)     if ($mappings_5);
printf("3 : mappings    : %7d\n",  $mappings_3) if ($mappings_on);
printf("3 :avg length   : %7.2f\n", $summed_3_length/ $mappings_3)   if ($mappings_3);
printf("3 :avg coverage : %7.2f\n", $summed_3_mean/ $mappings_3)     if ($mappings_3);

print "no mappings in $zero_mappings target region(s) \n";

__END__
#$off_target_profile->dump_profile();

print "--------------------------------------\n";
$profile->dump_profile();
print "--------------------------------------\n";

print "depth mean: ". $profile->mean_coverage() . "\n";
print "depth modal: ". $profile->modal_coverage() . "\n";


__END__ 

my $timer = EASIH::Timer->New();

my $stack = EASIH::Stacker->New;

my $START = 0;
my $END   = 1;
my $SUB   = 2;
my $LEVEL = 2;

my $VERBOSE_LEVEL = 0;

my @test_data = ([100,150], [98,101], [87,89], [175,180], [155, 170], [200,250], [185,190]);

@test_data = ([100,150], [200,250], [300, 350], [180,181], [165, 170],  [285,290], [155,156]);
@test_data = ([100,150], [50,95], [90,110], [90, 105]  );

my @stack;# = shift @test_data;

my $dprint = 2;

$timer->start();
my ($identical, $contained, $overlapped, $pushed) = (0,0,0, 0);

print "-"x60 . "\n";
print " "x25 . "START RUN\n";
print "-"x60 . "\n";

if (1) {

  my $features;
  while(<>) {
    chomp;
    
    next if (/^\z/);
    push @$features, [split("\t", $_)];
  }
  $stack->inserts_sorted($features, 1);
}
else {
  b_inserts_sorted(\@stack, \@test_data, 1);
}

$timer->stop();

my $fragment_cutoff = 0;

print "-"x60 . "\n";
print " "x25 . "ENDED RUN\n";
print " "x20 . "RUNTIME :: ".( $timer->report('s'))."\n";
print "-"x60 . "\n";

#print Dumper($stack->stack2array($fragment_cutoff));

$stack->print_stack();
my $profile =  $stack->profile();
print "\n";
$profile->dump_profile();
print "--------------------------------------\n";
#$profile->lower_resolution(20);
$profile->dump_profile();
print " Storing profile\n";
$profile->export_data('dil.1.profile');
#exit;

#my @outlines = $stack->stack2outline();
#map { print "$$_[ 0 ] --> $$_[ 1 ] |$$_[ 2 ]|\n"}  @outlines;
#print $stack->stack2bedGraph('chr1');
print "depth mean: ". $profile->mean_coverage() . "\n";
print "depth modal: ". $profile->modal_coverage() . "\n";
print "depth median: ". $profile->median_coverage() . "\n";
print "depth weighted mean: ". $profile->weighted_mean_coverage() . "\n";

#my @blocks = $stack->covered_regions();
#map{ print "--> @$_ \n"} @blocks;

#open (my $out_dump, "> stack_dump") || die "Could not open file 'stack_dump': $1\n";
#print $out_dump $stack->export_data('tyt');
#close $out_dump;

#print Dumper( $dd );
#print Dumper( Storable::thaw($dd));

#my $new_stack = EASIH::Stacker->New;

#$new_stack->import_data( $dd );
#print $new_stack->dump_stack() . "\n";
#print "MAX DEPTH: " . $stack->max_depth( 0) . "\n";

#print Dumper( flatten_stack( \@stack, $fragment_cutoff));

my($nodes, $max_depth, $troughs) = $stack->stats($fragment_cutoff );

print "Stack stats :: nodes= $nodes, max depth= $max_depth, troughs= $troughs \n";
