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
my $logic_name  = "illumina";
my $write       = 1;
my $verbose     = 0;
my $stacked     = 0;
my $min_score   = 5;
my $tab_output  = 0;

&GetOptions('dbhost:s'      => \$dbhost,
            'dbport:n'      => \$dbport,
            'dbname:s'      => \$dbname,
            'dbuser:s'      => \$dbuser,
            'dbpass:s'      => \$dbpass,
            'logic_name:s'  => \$logic_name,
            'analysis:s'    => \$logic_name,
	    'stacked'       => \$stacked,
	    'score:s'       => \$min_score,
            'write'         => \$write,
            'verbose'       => \$verbose,
	    'tab_output'    => \$tab_output,
           );

if (! $dbhost || !$dbname || ! $dbuser || !$dbpass || ! $logic_name ) {
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
#  print Dumper( $analysis );
$aa->store( $analysis);
}

my $chromosome_name = '';
my $mappings = 0;

my %seq_region_id_cache = ();

while (<STDIN>) { 

    
  chomp;
  
  my($q_name, $flags, $t_name, $start, 
     $mapq, $cigar, $mate_ref, 
     $mate_pos, $insert_size, $sequence,
     $quality, $opt_flags) = split("\t");

##  print "  next if ( $min_score && $min_score > $mapq ); \n";
  next if ( $min_score && $min_score > $mapq );
  
  my $seq_region_id = $seq_region_id_cache{ $t_name };
  if ( ! $seq_region_id ) {
    my $slice = $sa->fetch_by_region( 'chromosome', $t_name );
    $seq_region_id = $slice->get_seq_region_id;
    $seq_region_id_cache{ $t_name } = $seq_region_id;
  }
  my $end = $start + length($sequence);

  if ( $tab_output ) {

    print join("\t", '\N',$seq_region_id,$start, $end, 1, 1, $end - $start + 1, 1, $q_name, $analysis->dbID, $mapq, '\N','\N', $cigar, '\N', '\N', '\N', '\N') . "\n";

  }
  else {
    $daf_insert->execute($seq_region_id, $start, $end , 1, 
			 1, $end - $start + 1, 1, $q_name, $analysis->dbID, $mapq ) if ( $write);
    $daf_insert->finish();
    $mappings++;
  }
    
}
