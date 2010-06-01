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
use lib '/home/kb468/projects/e57/bioperl-live/';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::AssemblyMapper;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;


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

my $asma = $db->get_AssemblyMapperAdaptor();
my $csa  = $db->get_CoordSystemAdaptor();
my $sa   = $db->get_SliceAdaptor();

my $from_cs = $csa->fetch_by_name('chromosome', $from );
die "Unknown coord system: $from\n" if ( !$from_cs );
my $to_cs   = $csa->fetch_by_name('chromosome', $to);
die "Unknown coord system: $to\n" if ( !$to_cs );


my $mapper  = $asma->fetch_by_CoordSystems( $from_cs, $to_cs );


print "Input tab seperated data (chr,start,end) or regions (chr:start-end)\n";

while(<>) {
  chomp;
  my ( $chr, $start, $end ) = split(/\s+/);

  if ( ( ! $chr || ! $start || !$end ) && 
       $_ =~ /(.*?):(\d+)-(\d+)/) {

    ( $chr, $start, $end ) = ( $1, $2, $3);
  }
  elsif ( ! $chr || ! $start || !$end ) {
    print "cannot extract chr start end (tab seperated) from this line '$_'\n";
    next;
  }

  $chr =~ s/chr//;
  my @res = $mapper->map($chr, $start, $end, 1, $from_cs);
  foreach my $res ( @res ) {
    if ( $res->isa( 'Bio::EnsEMBL::Mapper::Coordinate' )) {
      my $chr_slice = $sa->fetch_by_seq_region_id($res->id);
      print "$chr, $start, $end --> " . join("\t", $chr_slice->seq_region_name, $res->start, $res->end ) . "\n";
    }
  }
}



