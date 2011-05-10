package EASIH::Toolbox;
# 
# Misc toolbox functions that does not fit elsewhere
# 
# 
# Kim Brugger (18 Jan 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Storable;
use File::Temp;
use Time::HiRes;
use Carp;

# 
# 
# 
# Kim Brugger (20 Sep 2010)
sub version {


  my $libdir = $INC{ 'EASIH/Toolbox.pm'};

  my $VERSION   = "unknown";

  if ($libdir && $libdir =~ /(.*)\//) {
    $VERSION = `cd $1; git describe --always --dirty 2> /dev/null`;
  }
  else {
    $VERSION = `git describe --always --dirty 2> /dev/null`;
  }
  $VERSION ||= "Unknown";

  chomp( $VERSION );

  return $VERSION;
}

1;




