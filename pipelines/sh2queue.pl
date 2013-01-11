#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (10 Oct 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;


my $opts = 'N:l:';
my %opts;
getopts($opts, \%opts);

my $qsub_params = "-cwd ";
foreach my $k ( keys %opts ) {
  $qsub_params .= " -$k $opts{ $k } ";
}


while(<>) {
  system "echo \" $_ \" | qsub $qsub_params \n"
}
