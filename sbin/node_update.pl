#!/usr/bin/perl 
# 
# script for updating the nodes. Should probably have done this in bash...
# 
# 
# Kim Brugger (25 Jul 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


my $selection_file = "/scratch2/systems/selections";
exit if ( ! -e $selection_file);
system "cat $selection_file | dpkg --set-selections ";

# always hold the kernel updates, just do it here to make sure that nothing slips through
system "dpkg --get-selections | grep linux | egrep '(firmware|headers|image|server)' | perl -pe 's/install/hold/' | sudo dpkg --set-selections ";

system "apt-get -y dselect-upgrade "



