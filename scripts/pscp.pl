#!/usr/bin/perl 
# 
# parallel scp, for speed!
# 
# 
# Kim Brugger (24 Nov 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my $debug = 0;
$debug = 1;

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

use EASIH::Parallel;

my %opts;
getopts("n:s:h", \%opts);


my (@files, $dest);

$dest = pop(@ARGV);
@files = @ARGV;


foreach my $file ( @files ) {

  EASIH::Parallel::job_push(\&scp, $file, $dest);
}

#EASIH::Parallel::dump_jobs();
#exit;
EASIH::Parallel::max_nodes( $opts{n}) if ( $opts{n} );
EASIH::Parallel::run_parallel( );



# 
# 
# 
# Kim Brugger (24 Nov 2011)
sub scp {
  my ($file, $dest) = @_;
  system "scp $file $dest\n";
}

