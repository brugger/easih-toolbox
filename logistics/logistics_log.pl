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
use EASIH::Logistics;

my $runfolder = shift || "";



if ( !$runfolder ) {
  my @data = EASIH::Logistics::fetch_offloaded_files();
#  print Dumper( \@data );
  foreach my $line ( @data ) {
    print "$$line[0] --> $$line[1]\n";
  }
}
else {
  my @data = EASIH::Logistics::fetch_files_from_rundir( $runfolder );
  foreach my $line ( @data ) {
#    print "$runfolder --> $line\n";
  }

  my $printstring = EASIH::Logistics::runfolder_log( $runfolder );

  print $printstring;
}
