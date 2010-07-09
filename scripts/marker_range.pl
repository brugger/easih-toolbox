#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (18 Jun 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/usr/local/lib/ensembl/modules/';
use lib '/usr/local/lib/bioperl/';
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor;
use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MarkerFeature;
use Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor;

use Getopt::Std;

my %opts;
getopts('s:h', \%opts);
usage() if ( $opts{h});

my $species     = $opts{s} || "human";
my $host        = 'ensembldb.ensembl.org';
my $user        = 'anonymous';

# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);
# get variation adaptors
my $ma  = $reg->get_adaptor($species, 'core', 'Marker');
my $mfa = $reg->get_adaptor($species, 'core', 'MarkerFeature');

while( $_ = shift) {
  chomp;

  if ( ! /\w:\w/ ) {
    print "$_ is not in the expected marker:marker format\n";
    next;
  }

  my ( $first, $second) = split(':');
  my ($first_chr, $first_start, $first_end)    = fetch_marker_postion( $first );
  my ($second_chr, $second_start, $second_end) = fetch_marker_postion( $second );

  die "The markers are on two different chromosomes ($first_chr, $second_chr) \n" if ( $first_chr ne $second_chr);

  print "$first:$second ==> ";

  if ( $first_start > $second_start) {
    print "$first_chr:$second_start-$first_end";
  }
  else {
    print "$first_chr:$first_start-$second_end";
  }
  print "\n";

}



# 
# 
# 
# Kim Brugger (18 Jun 2010)
sub fetch_marker_postion {
  my ( $name ) = shift;


  my @markers = @{$ma->fetch_all_by_synonym($name)};
  die "More than one synonymr $name !!\n" if ( @markers > 1 );
  
  my @features = @{$mfa->fetch_all_by_Marker($markers[0])};
  die "More than one feature for $name !!\n" if ( @features > 1 );
  
  return ($features[0]->seq_region_name(),$features[0]->seq_region_start(), $features[0]->seq_region_end());
}




# 
# 
# 
# Kim Brugger (09 Jul 2010)
sub usage {

  print "marker_range.pl finds the chromosome positions between two markers. The script uses the current ensembl version of the species of interest.\n";
  print "Usage: marker_range.pl <options> markers. The markers can be given on the command line seperated by ':'.\n";
  exit;
}


=pod

=head1 SYNOPSIS

marker_range.pl finds the chromosome positions between two markers. The script uses the current ensembl version of the species of interest.

=head1 OPTIONS

Usage: marker_range.pl <options> markers. The markers can be either given on the command line seperated by ':'.

=cut

