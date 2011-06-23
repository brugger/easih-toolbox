#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (23 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


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
use EASIH::QC;
use EASIH::QC::db;


my $sid = EASIH::QC::db::add_sample('Z980001', 'Z98');
my $pid = EASIH::QC::db::fetch_project_id('Z98');
print "$pid --> $sid\n";

