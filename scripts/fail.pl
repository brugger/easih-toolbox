#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (10 Aug 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $sleep = shift || int(rand(100)+10);

print "Sleeping $sleep seconds\n";

sleep( $sleep );

exit 1;
