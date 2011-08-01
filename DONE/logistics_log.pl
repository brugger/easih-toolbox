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

my $runfolder = shift || "";


if ( !$runfolder ) {
  my @data = EASIH::DONE::fetch_runs();
#  print Dumper( \@data );
  foreach my $line ( @data ) {
    print "$$line[1] --> $$line[2]\n";
  }
}
else {
  my @data = EASIH::DONE::fetch_files_from_run( $runfolder );
  foreach my $line ( @data ) {
  }

  my $printstring = EASIH::DONE::run_log( $runfolder );

  print $printstring;
}
