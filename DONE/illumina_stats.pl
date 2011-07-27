#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
my $DYNAMIC_LIB_PATHS = 1;
BEGIN {
  if ( $DYNAMIC_LIB_PATHS ) {
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    use lib '/home/kb468/easih-toolbox/modules/';
  }

}

use EASIH;
use EASIH::DONE;

my $run_folder = shift;

if ( ! $run_folder) {
  my @runs = EASIH::DONE::fetch_runs();
  foreach my $run ( @runs ) {
    print join("\t", $$run[1], $$run[2], "\n");
  }
}
else {
  my $rid = EASIH::DONE::fetch_run_id($run_folder);

  my @lanes  = EASIH::DONE::fetch_illumina_lane_stats_by_rid($rid);
  my @mplexs = EASIH::DONE::fetch_illumina_multiplex_stats_by_rid($rid);

  print join("\t", "lane", "Sample", "total", "PF", "Stats", "\n");
  print "="x50 ."\n";
  foreach my $lane ( @lanes) {
    print join("\t", "lane $$lane[2].$$lane[3]", $$lane[4], $$lane[5], $$lane[6], sprintf("%.2f %%", $$lane[6]*100/$$lane[5]), "\n" );
    if ( $$lane[4] eq "MULTIPLEXED") {
      foreach my $mplex (@mplexs ) {
	if ( $$lane[2] == $$mplex[2] ) {
	  print join("\t", $$mplex[3], $$mplex[4], $$mplex[6], "$$mplex[5] %", "\n");
	}
      }
    }

  }

}
