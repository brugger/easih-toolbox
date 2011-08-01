#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (15 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


use strict;
use warnings;
use Data::Dumper;


# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 1;
  if ( $DYNAMIC_LIB_PATHS ) {
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules" && ! grep /^$path\/modules/, @INC);
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    push @INC, '/home/kb468/easih-toolbox/modules/';
  }
}


use EASIH;
use EASIH::DONE;

my $runfolder = shift || usage();
my $rid = EASIH::DONE::fetch_run_id( $runfolder );

my $last_status = EASIH::DONE::fetch_latest_offloading_status( $rid );

die "no run with the name of '$runfolder'\n" if (! $last_status );

EASIH::DONE::add_offloading_status($rid, 
					"RETRY_OFFLOAD");




# 
# 
# 
# Kim Brugger (22 Jun 2011)
sub usage {
  print "usage: $0 <run-folder/name>\n";
  exit -1;
}
