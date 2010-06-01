#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (27 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


use lib '/home/kb468/projects/e57/ensembl/modules/';
use lib '/home/kb468/projects/e57/ensembl-variation/modules/';
use lib '/home/kb468/projects/e57/bioperl-live/';

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);


my $variation_adaptor = $registry->get_adaptor(
	'human',	# species
	'variation',	# database
	'variation'	# object type
);

my $sa = $registry->get_adaptor( 'Human', 'Core', 'Slice' );


my $chr   = 13;
my $start = 32927286;
my $end   = 32930082;
# $end   = 32927286;
$start = 32927607;
$end = $start;

my $slice = $sa->fetch_by_region('chromosome', $chr, $start, $end, 1);

my $genes = $slice->get_all_Genes();
foreach my $gene ( @$genes ) {

 print("gene info :: "
      . join( ":", map { $gene->$_ } qw(stable_id external_name start end description ) )
      . "\n" );


}

my $transcripts = $slice->get_all_Transcripts();
foreach my $transcript ( @$transcripts ) {

}

my $vars = $slice->get_all_VariationFeatures();
foreach my $var ( @$vars ) {

  print("variation info :: "
      . join( ":", map { $var->$_ } qw(variation_name allele_string start end ) )
      . "\n" );


}


