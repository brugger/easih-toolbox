#!/usr/bin/perl 
# 
# Transform human GRCh37/hg18 coordinates into human NCBI36/hg18 coordinates 
# 
# 
# Kim Brugger (27 May 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


use lib '/software/lib/ensembl-variation/modules/';
use lib '/software/lib/ensembl-functgenomics/modules/';
use lib '/software/lib/ensembl/modules/';
use lib '/software/lib/bioperl/';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::AssemblyMapper;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;


my $species     = "human";
my $to          = 'NCBI36';
my $from        = 'GRCh37';
my $host        = 'ensembldb.ensembl.org';
my $user        = 'anonymous';

# get registry
my $reg = 'Bio::EnsEMBL::Registry';
#$reg->load_registry_from_db(-host => $host,-user => $user);
$reg->load_registry_from_db(-host => "mgpc17",-user => "easih_ro", -NO_CACHE => 0);

my $asma = $reg->get_adaptor($species, 'core', 'AssemblyMapper');
my $csa  = $reg->get_adaptor($species, 'core', 'CoordSystem');
my $sa   = $reg->get_adaptor($species, 'core', 'Slice');


my $from_cs = $csa->fetch_by_name('chromosome', $from );
die "Unknown coord system: $from\n" if ( !$from_cs );
my $to_cs   = $csa->fetch_by_name('chromosome', $to);
die "Unknown coord system: $to\n" if ( !$to_cs );


my $mapper  = $asma->fetch_by_CoordSystems( $from_cs, $to_cs );

my $vcf = 0;

#print "Input tab seperated data (chr,start,end) or regions (chr:start-end)\n";

while(<>) {

  if (/#/) {
    print;
    next;
  }

  chomp;
  my ( $chr, $start, $end, $rest ) = split(/\s+/);

  if ( ( ! $chr || ! $start || !$end ) && 
       $_ =~ /(.*?):(\d+)-(\d+)/) {
    ( $chr, $start, $end ) = ( $1, $2, $3);
  }
  # This is a vcf format file..
  elsif ( $end && $end !~ /^\d+\z/) {
    $vcf = 1;
    $end = $start;
  }
  elsif ( ! $chr || ! $start || !$end ) {
    print "cannot extract chr start end (tab seperated) from this line '$_'\n";
    next;
  }

  $chr =~ s/chr//;
  my @res;

  eval { @res = $mapper->map($chr, $start, $end, 1, $from_cs)};
  foreach my $res ( @res ) {
    if ( $res->isa( 'Bio::EnsEMBL::Mapper::Coordinate' )) {
      my $chr_slice = $sa->fetch_by_seq_region_id($res->id);
#      print "$chr, $start, $end --> " . join("\t", $chr_slice->seq_region_name, $res->start, $res->end ) . "\n";
      if ( $vcf ) {
	my @f = split("\t");
	($f[0], $f[1]) = ($chr, $start);
	print join("\t", @f) . "\n";
      }
      else {
	print join("\t", $chr_slice->seq_region_name, $res->start, $res->end, $rest ) . "\n";
      }
    }
    elsif ( $res->isa( 'Bio::EnsEMBL::Mapper::Gap' )) {
      print STDERR "$chr:$start-$end is in a GAP\n";
    }
    else {
      print STDERR "Could not transform $chr:$start-$end\n";
    }
  }
}



