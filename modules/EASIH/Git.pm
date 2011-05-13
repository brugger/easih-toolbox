package EASIH::Git;

# 
# Some simple git interaction scripts.
# 
# If this file is included and it is not in the master branch a warning will come.
# 
# Kim Brugger (13 May 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;

# 
# Find the git version of the program using the EASIH::Git package otherwise the package
# given as argument 
#
# Kim Brugger (20 Sep 2010)
sub branch {
  my ( $module ) = @_;

  my $file = $0;

  if ( $module ) {
    $module =~ s/\//::/;
    $module =~ s/\;//;
    $module .= ".pm";
    $file = $INC{ $module} 
  }

  my $branch   = "unknown";

  if ($file && $file =~ /(.*)\//) {
    $branch = `cd $1; git branch | egrep ^\* 2> /dev/null`;
  }
  else {
    $branch = `git branch | egrep ^\* 2> /dev/null`;
  }
  $branch ||= "Unknown";

  chomp( $branch );
  $branch =~ s/^\*//;

  return $branch;
}


# 
# Find the git version of the program using the EASIH::Git package otherwise the package
# given as argument 
#
# Kim Brugger (20 Sep 2010)
sub version {
  my ( $module ) = @_;

  my $file = $0;

  if ( $module ) {
    my $mfile = $module;
    $mfile =~ s/::/\//;
    $mfile =~ s/\;//;
    $mfile .= ".pm";
    $file = $INC{ $mfile };
    
    return "Cannot find version for '$module' it is not loaded\n" if ( ! $file );

  }

  print STDERR "$file || $module\n";

  my $VERSION   = "unknown";

  if ($file && $file =~ /(.*)\//) {
    $VERSION = `cd $1; git describe --always --dirty 2> /dev/null`;
  }
  else {
    $VERSION = `git describe --always --dirty 2> /dev/null`;
  }
  $VERSION ||= "Unknown";

  chomp( $VERSION );

  return $VERSION;
}



BEGIN { 
  
  my $branch = branch();
  if ( $branch ne "master" ) {

    print STDERR "="x65 . "\n";
    print STDERR "="x5 ." THIS IS NOT THE MASTER BRANCH, NOT FOR PRODUCTION USE " . "="x5 . "\n";
    print STDERR "="x65 . "\n\n";

  }

  
  
}

1;



