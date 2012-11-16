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

  my $run_config = XMLin($in);

#  print Dumper( $config );

  my %config = ( raw => $run_config);
#  my %config = ( raw => '');

  

  $config{ is_indexed } = 0;
  $config{ is_indexed } = 1 if ($$run_config{ 'Run' }{ 'RunParameters' }{ 'Barcode'} );


  if (ref($$run_config{ 'Run' }{'TileSelection'}{ 'Lane' }) eq "ARRAY") {
    map { push @{$config{lanes}}, $$_{ 'Index' } } @{$$run_config{ 'Run' }{ 'TileSelection' }{ 'Lane' }};
  }
  else {
    push @{$config{lanes}}, $$run_config{ 'Run' }{ 'TileSelection' }{ 'Lane' }{ 'Index' };
  }

#  print Dumper($$run_config{ 'Run' });

  if (ref($$run_config{ 'Run' }{ 'RunParameters' }{ 'Reads' }) eq "ARRAY") {
    map { push @{$config{reads}}, $$_{ 'Index' } } @{$$run_config{ 'Run' }{ 'RunParameters' }{ 'Reads' }};
  }
  else {
    push @{$config{reads}}, $$run_config{ 'Run' }{ 'RunParameters' }{ 'Reads' }{ 'Index' };
  }
    

#  die Dumper ( \%config );
  
    

  return \%config;

#  return \%res;
}


1;



