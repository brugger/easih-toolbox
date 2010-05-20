#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (17 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $sleep_time = shift || 300;

srand;
$sleep_time = int(rand( $sleep_time ));

print "Starting sleep for $sleep_time seconds\n";
sleep($sleep_time);
print "Done sleeping.\n";

