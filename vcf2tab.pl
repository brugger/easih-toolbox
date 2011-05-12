#!/usr/bin/perl 
# 
# reads in the snp information into the database
# 
# 
# Kim Brugger (26 Nov 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

require DBI;
my $dbi = DBI->connect('DBI:mysql:dbsnp132:mgpc17', 'easih_ro') || die "Could not connect to database: DBI::errstr";

my $hgmd = shift || die " $0 hgmd-file dbsnp-VCF-file\n";
my $vcf  = shift || die "needs a dbsnp vcf file\n";

my $ref_id_GRCh37 = check_ref( "GRCh37" );
my $ref_id_hg18   = check_ref( "hg18" );



use lib '/software/lib/e62/ensembl-variation/modules/';
use lib '/software/lib/e62/ensembl-functgenomics/modules/';
use lib '/software/lib/e62/ensembl/modules/';
use lib '/software/lib/bioperl/';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::AssemblyMapper;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;


my $species     = "human";
my $from        = 'NCBI36';
my $to          = 'GRCh37';

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => "mgpc17",-user => "easih_ro", -NO_CACHE => 0);

my $asma = $reg->get_adaptor($species, 'core', 'AssemblyMapper');
my $csa  = $reg->get_adaptor($species, 'core', 'CoordSystem');
my $sa   = $reg->get_adaptor($species, 'core', 'Slice');

my $from_cs = $csa->fetch_by_name('chromosome', $from );
die "Unknown coord system: $from\n" if ( !$from_cs );
my $to_cs   = $csa->fetch_by_name('chromosome', $to);
die "Unknown coord system: $to\n" if ( !$to_cs );
my $mapper  = $asma->fetch_by_CoordSystems( $from_cs, $to_cs );


my $hgmds = undef;#readin_hgmd($hgmd);

open (my $flags, " > flags.txt") || die "Could not open 'flags.txt': $!\n";
open (my $snps,  " > snp.txt"  ) || die "Could not open 'flags.txt': $!\n";
open (my $meta,  " > meta.txt" ) || die "Could not open 'flags.txt': $!\n";

open( my $v, $vcf) || die "Could not open '$vcf': $!\n";
while (<$v>) {

  chomp;

  if ( /##INFO=\<ID=(\w+?),.*Description=\"(.*?)\"/ ) {
    my ($id, $full) = ($1, $2);
    print $flags join("\t", $id, $full, "\n");
  }
  next if (/#/);

  chomp;
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
  my ($chr, $pos, $rs, $ref_base, $alt_base, undef, undef, $flags, $format) = split("\t");
  
  next if ( $chr eq "PAR");

  $chr =~ s/chr//;

#  exit;

  my $class = "UNKNOWN";
  $flags =~ /VC=(.*?)\;/;
  $class = $1;

  my @flags = split(";", $flags);
  $flags = join(";", grep (/AF=|CDA=|CLN=|dbSNPBuildID=|G5=|G5A=|HD=|KGPilot|KGPROD=|KGVAL=|MUT=|NS=|OM=|PM=|VLD/, @flags));

  
  my $HGMD = "N";
  $HGMD = "Y" if ( $$hgmds{$rs} || $$hgmds{$chr}{$pos});

  my @rss = split(";", $rs);
  foreach ( @rss ) {
    my @f = ($_, $chr, $pos, $ref_base, $alt_base, $class, $HGMD, $flags);
    print $snps join("\t", @f) . "\n";
  }

}
close($v);


# 
# 
# 
# Kim Brugger (03 Dec 2010)
sub remap {
  my ($chr, $start, $end) = @_;
  
  $end ||= $start;

  $chr =~ s/chr//;
  my @res = $mapper->map($chr, $start, $end, 1, $from_cs);
  foreach my $res ( @res ) {
    if ( $res->isa( 'Bio::EnsEMBL::Mapper::Coordinate' )) {
      my $chr_slice = $sa->fetch_by_seq_region_id($res->id);
      return($chr_slice->seq_region_name, $res->start, $res->end );
    }
    elsif ( $res->isa( 'Bio::EnsEMBL::Mapper::Gap' )) {
    }
    else {
    }
  }
  return (undef, undef, undef);
}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub readin_hgmd {
  my ( $file ) = @_;
  
  my %res;
  my %done;
  open ( my $f, $file) || die "Could not open '$file': $!\n";
  while (<$f>) {

    next if (/#/ || /^\/\// || /^Type/);
    
    chomp;
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    my @f = split("\t");
    my $rs = $f[3];
    my ($chr, $pos) = split(":", $f[4]);
    $chr =~ s/chr//;
    $pos ||= -1;

    if ( $rs eq 'null' && $chr ne "null") {
      next if ($done{$ref_id_hg18}{$chr}{$pos});
      my @f = ('', $chr, $pos, '', '', $ref_id_hg18, "Y", '', '');
#      print join("\t", @f) . "\n";

      $done{$ref_id_hg18}{$chr}{$pos}++;

      ($chr, $pos) = remap($chr, $pos);
      next if ( ! $chr);
      $chr =~ s/chr//;
      @f = ('', $chr, $pos, '', '', $ref_id_GRCh37, "Y", '', '');
#      print join("\t", @f) . "\n";

      next;
    }

    next if ($done{$rs});
    $res{$rs} = "Y";
    $done{$rs}++;
    next if ($chr eq 'null');
    ($chr, $pos) = remap($chr, $pos);
    next if ( ! $chr);
    $chr =~ s/chr//;

    $res{$chr}{$pos} = "Y";
    
  }
  close ($f);

  return \%res;
}


# 
# 
# 
# Kim Brugger (26 Nov 2010)
sub check_ref {
  my ($name) = @_;

  my $sth = $dbi->prepare("select ref_id from reference where name='$name' or alias like '%$name%';");
  $sth->execute();
  my @entries =  $sth->fetchrow_array();
  if (@entries) {
    return $entries[0];
  }

  $sth = $dbi->prepare("INSERT INTO reference (name) VALUES ('$name');");
  $sth->execute();
  $sth = $dbi->prepare("select ref_id from reference where name='$name' or alias like '%$name%';");
  $sth->execute();
  @entries =  $sth->fetchrow_array();
  if (@entries) {
    return $entries[0];
  }

  die "Could not find or insert reference '$name'\n";
}
