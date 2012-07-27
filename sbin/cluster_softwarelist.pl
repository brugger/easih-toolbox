#!/usr/bin/perl 
# 
# dumps a list of installed files from the server, to be collected by the cluster nodes...
# 
# 
# Kim Brugger (25 Jul 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $selection_file = "/scratch2/systems/selections";
my $tmp_selection_file = "/scratch2/systems/selections.tmp";


my @selections = split("\n", `dpkg --get-selections`);

@selections = grep(! /(firmware|header|image|gridengine|linux-server|bind9|dhcp3-server)/, @selections);

open (my $o, "> $tmp_selection_file") || die "Could not write to '$tmp_selection_file': $!\n";
print $o join("\n", @selections) . "\n";
close $o;

if (! -e $selection_file )  {
  system "mv $tmp_selection_file $selection_file";
  exit;
}
   

my $changed_lines = `cat  $selection_file $tmp_selection_file | sort | uniq -u | wc -l`;
chomp $changed_lines;

if ( $changed_lines >= 1 ) {
  system "mv $tmp_selection_file $selection_file";
}  

