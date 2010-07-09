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


my %stacks;

my $chromosome_name = '';
my $mappings = 0;

while (<>) { 
  
  chomp;

  my ($seq_region, undef, undef, $start, $end, $strand, undef) = split("\t");
  my $slice = $sa->fetch_by_region( 'chromosome', $seq_region );

  my $seq_region_id;
  eval { $seq_region_id = $slice->get_seq_region_id};
  next if ( ! $seq_region_id );

  $daf_insert->execute($slice->get_seq_region_id, $start, $end , 1, 
		       1, $end - $start + 1, 1, "no_name", $analysis->dbID, $end - $start + 1) if ( $write);
  $daf_insert->finish();
  $mappings++;

}



print "Inserted $mappings mapped regions \n"; 


__END__





  my ($start, $end) = (0,0);
for ( ;$end < @chromosome; $end++) {
  if ( $chromosome[$end] && !$start) {
    $start = $end;
    next;
  }
  elsif (!$chromosome[$end] && $start) {
#    $daf_insert->execute($slice->dbID, $start, $end - 1 , 1, 1, $end - $start + 1, 1, "no_name", $analysis->dbID, $end - $start + 1) if ( $write);
#    $daf_insert->finish();

    print " $start -> $end \n";


    $mappings++;
    $start = undef;
  }
}

if ( $start ) { 
  $daf_insert->execute($slice->dbID, $start, $end - 1 , 1, 1, $end - $start + 1, 1, "no_name", $analysis->dbID, $end - $start + 1) if ( $write);
  $daf_insert->finish();
    $mappings++;
    print " $start -> $end \n";
}

