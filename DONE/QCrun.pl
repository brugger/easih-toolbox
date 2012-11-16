#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (18 Sep 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Getopt::Std;

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

my $easih_toolbox = '/software/installed/easih-toolbox/';

my %opts;
getopts('hr:', \%opts);

usage() if ( $opts{h} );

my $rid = EASIH::DONE::fetch_run_id( $opts{'r'} );

print STDERR "rid : $rid\n";


my %fqfiles = EASIH::DONE::fetch_files_from_rid_by_sample( $rid );

#print Dumper( \%fqfiles );
foreach my $sample (keys %fqfiles)  {
  my $cmd   = "$easih_toolbox/DONE/QC_report.pl -p ILLUMINA -r ";
  $cmd     .= " -1 $fqfiles{ $sample }[0] ";
  $cmd     .= " -2 $fqfiles{ $sample }[1]" if ($fqfiles{ $sample }[1]);

  print "$cmd \n";

}



# 
# 
# 
# Kim Brugger (18 Sep 2012)
sub usage {
  print "QCrun.pl : -r<un dir/id> -h[elp]\n";
  exit -1;


}

