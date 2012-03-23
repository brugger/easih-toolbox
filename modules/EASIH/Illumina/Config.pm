package EASIH::Illumina::Config;
# 
# For pulling data out of the BaseCalls/config.xml file
# 
# 
# Kim Brugger (17 Nov 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use XML::Simple;

# 
# 
# 
# Kim Brugger (17 Nov 2011)
sub readin {
  my ($in) = @_;
  
  $in .= '/config.xml' if ( -d $in );

  my $config = XMLin($in);

#  print Dumper( $config );

  return $config;

#  return \%res;
}


1;



