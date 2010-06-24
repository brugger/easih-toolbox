#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (24 Jun 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


use lib '/home/kb468/projects/easih-flow/modules';
use lib '/home/kb468/easih-pipeline/modules';
use EASIH::JMS;

my $freeze_file = shift || usage();


EASIH::JMS::restore_state($freeze_file);
print EASIH::JMS::job_report( 1 );
EASIH::JMS::run_stats();    
EASIH::JMS::memory_stats();    





# 
# 
# 
# Kim Brugger (24 Jun 2010)
sub usage {
  
  $0 =~ s/.*\///;
  die "USAGE: $0 pipeline-file\n";
  
}

