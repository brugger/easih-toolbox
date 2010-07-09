#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (26 Mar 2010), contact: kim.brugger@easih.ac.uk

use lib '/home/kb468/projects/e57/ensembl/modules/';
use lib '/home/kb468/projects/e57/ensembl-pipeline/modules/';
use lib '/home/kb468/projects/e57/bioperl-live/';


use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;


my $dbhost      = "localhost";
my $dbport      = 3306;
my $dbname      = "human_GRCh37";
my $dbpass      = "ensembl";
my $dbuser      = "ensadmin";
my $logic_name  = "454";
my $blat_file   = "454.blat.filtered";
my $write       = 1;
my $verbose     = 0;

&GetOptions('dbhost:s'      => \$dbhost,
            'dbport:n'      => \$dbport,
            'dbname:s'      => \$dbname,
            'dbuser:s'      => \$dbuser,
            'dbpass:s'      => \$dbpass,
            'logic_name:s'  => \$logic_name,
            'analysis:s'    => \$logic_name,
            'write'         => \$write,
            'verbose'       => \$verbose,
            'blat_file:s'   => \$blat_file,
            'blat:s'        => \$blat_file,
           );

if (! $dbhost || !$dbname || ! $dbuser || !$dbpass || ! $logic_name || !$blat_file) {
  system "/usr/bin/perldoc $0";
  exit;
}


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                                       -dbname => $dbname,
                                                       -host   => $dbhost,
                                                       -user   => $dbuser,
                                                       -port   => $dbport,
                                                       -pass   => $dbpass,
                                                       );

my $dbc  = $db->dbc;
my $sa   = $db->get_SliceAdaptor();
my $dafa = $db->get_DnaAlignFeatureAdaptor();
my $aa   = $db->get_AnalysisAdaptor();

my $daf_insert = $dbc->prepare( qq{ insert into dna_align_feature (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_start, hit_end, hit_strand, hit_name, analysis_id, score) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)});

my $analysis = $aa->fetch_by_logic_name( $logic_name );
if (! $analysis ) {
  $analysis = new Bio::EnsEMBL::Analysis( -logic_name => $logic_name );
  print Dumper( $analysis );
  $aa->store( $analysis);
}


open (my $blat, $blat_file) || die "Could not open '$blat_file': $!\n";

while (<$blat>) { 
  
  chomp;

  my($match, $mis_match, $rep, $Ns, 
     $q_ins, $q_ins_bases, $t_ins, 
     $t_ins_bases, $strand, $q_name, 
     $q_length, $q_align_start, $q_align_end, 
     $t_name, $t_size, $t_align_start, 
     $t_align_end, $blocks, $block_lengths, 
     $q_block_starts, $t_block_starts) = split("\t");


  my $slice = $sa->fetch_by_region( 'chromosome', $t_name );

  $daf_insert->execute($slice->dbID, $t_align_start, $t_align_end, $strand, $q_align_start, $q_align_end, 1, $q_name, $analysis->dbID, $match) if ( $write);
  $daf_insert->finish();
}




