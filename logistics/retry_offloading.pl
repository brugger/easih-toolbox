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

use DBI;


# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment. Needs to be prior to the use of EASIH* modules.
BEGIN {
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

use EASIH;
use EASIH::Logistics;

my $runfolder = shift || usage();

my $last_status = EASIH::Logistics::fetch_latest_run_folder_status($runfolder);

die "no run with the name of '$runfolder'\n" if (! $last_status );

my $platform = EASIH::Logistics::fetch_platform($runfolder);


EASIH::Logistics::add_run_folder_status($runfolder, 
					$platform, 
					"RETRY_OFFLOAD");




# 
# 
# 
# Kim Brugger (22 Jun 2011)
sub usage {
  print "usage: $0 <run-folder/name>\n";
  exit -1;
}
