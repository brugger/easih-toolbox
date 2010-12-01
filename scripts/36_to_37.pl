#!/usr/bin/perl 
# 
# Transform human NCBI36/hg18 coordinates into human GRCh37/hg18 coordinates 
# 
# 
# Kim Brugger (27 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

use lib '/usr/local/lib/ensembl-variation/modules/';
use lib '/usr/local/lib/ensembl-functgenomics/modules/';
use lib '/usr/local/lib/ensembl/modules/';
use lib '/usr/local/lib/bioperl/';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::AssemblyMapper;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

my $species     = "human";
my $from        = 'NCBI36';
my $to          = 'GRCh37';
my $host        = 'ensembldb.ensembl.org';
my $user        = 'anonymous';

# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);

my $asma = $reg->get_adaptor($species, 'core', 'AssemblyMapper');
my $csa  = $reg->get_adaptor($species, 'core', 'CoordSystem');
my $sa   = $reg->get_adaptor($species, 'core', 'Slice');

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
#      print "$chr, $start, $end --> " . join("\t", $chr_slice->seq_region_name, $res->start, $res->end ) . "\n";
      print join("\t", $chr_slice->seq_region_name, $res->start, $res->end ) . "\n";
    }
    elsif ( $res->isa( 'Bio::EnsEMBL::Mapper::Gap' )) {
      print STDERR "$chr:$start-$end is in a GAP\n";
    }
    else {
      print STDERR "Could not transform $chr:$start-$end\n";
    }
  }
}



