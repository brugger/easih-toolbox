#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (17 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $sleep_time = shift || 10;
my $fail_rate = 50;

srand;
$sleep_time = int(rand( $sleep_time ));

#print "Starting sleep for $sleep_time seconds\n";
sleep($sleep_time);

if ( $fail_rate > int(rand(100))) {
  print "------- Failed -----------\n";
  exit 1;
}
#print "Done sleeping.\n";

