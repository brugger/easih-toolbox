#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (27 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use Getopt::Long;

use lib '/home/kb468/projects/e57/ensembl/modules/';
use lib '/home/kb468/projects/e57/ensembl-variation/modules/';
use lib '/home/kb468/projects/e57/bioperl-live/';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::AssemblyMapper;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;

use Bio::EnsEMBL::Variation::Variation;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;

my $dbhost      = "localhost";
my $dbport      = 3306;
my $dbname      = "human_e57";
my $dbpass      = "";
my $dbuser      = "easih_ro";
my $from        = 'NCBI36';
my $to          = 'GRCh37';

&GetOptions('dbhost:s'      => \$dbhost,
            'dbport:n'      => \$dbport,
            'dbname:s'      => \$dbname,
            'dbuser:s'      => \$dbuser,
	    'dbpass:s'      => \$dbpass,
	    'from:s'        => \$from,
	    'to:s'          => \$to,
           );



my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $dbname, -host   => $dbhost,
                                                       -user   => $dbuser, -port   => $dbport,
                                                       -pass   => $dbpass,                    );

my $sa   = $db->get_SliceAdaptor();


my $chr   = 13;
my $start = 32927286;
my $end   = 32930082;

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
  
  print   print $var->source(), ':',$var->name(), ".",$var->version,"\n";


}


