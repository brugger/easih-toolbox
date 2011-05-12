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
my $dbi = DBI->connect('DBI:mysql:dbsnp_132', 'easih_admin', 'easih') || die "Could not connect to database: DBI::errstr";


my $ref_id_GRCh37 = check_ref( "GRCh37" );
my $ref_id_hg18   = check_ref( "hg18" );

my %populations;

my $sth = $dbi->prepare("select * from populations;");
$sth->execute();
while (my $res = $sth->fetchrow_hashref()) {
  $populations{ $res->{short}} = $res->{pop_id};
}

my $pop;

while (<>) {

  if (/##dbSNP_LOC_POP_ID=HapMap-(\w+)/i) {
    $pop = $1;
    next;
  }

  next if (/#/);

  chomp;
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
  my (undef, undef, $rs, undef, undef, undef, undef, $flags, undef) = split("\t");

  my $allele_freq = $1 if ($flags  =~ /AF=(\d+\.\d+)/);
  my $sample_size = $1 if ($flags  =~ /NS=(\d+)/);
  
  die  "could not find pop flag\n" if ( !$pop );
  check_or_add_pop( $pop );

  if ( ! $sample_size || !$allele_freq || ! $populations{ $pop } || ! $rs) {
    print STDERR "$_\n";
    print STDERR "$rs, $populations{ $pop }, $sample_size, $allele_freq ($pop)\n";
    die;
  }
  
  print join("\t", $rs, $populations{ $pop }, $sample_size, $allele_freq) . "\n";

}

# 
# 
# 
# Kim Brugger (26 Nov 2010)
sub check_ref {
  my ($name) = @_;

  my $sth = $dbi->prepare("select ref_id from reference where name='$name';");
  $sth->execute();
  my @entries =  $sth->fetchrow_array();
  if (@entries) {
    return $entries[0];
  }

  $sth = $dbi->prepare("INSERT INTO reference (name) VALUES ('$name');");
  $sth->execute();
  $sth = $dbi->prepare("select ref_id from reference where name='$name';");
  $sth->execute();
  @entries =  $sth->fetchrow_array();
  if (@entries) {
    return $entries[0];
  }

  die "Could not find or insert reference '$name'\n";
}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub check_or_add_pop {
  my ($name) = @_;

  return if ($populations{ $name });

  $sth = $dbi->prepare("INSERT INTO populations (short, descr) VALUES ('$name', 'undefined');");
  $sth->execute();
  $sth = $dbi->prepare("select pop_id from populations where short='$name';");
  $sth->execute();
  my @entries =  $sth->fetchrow_array();
  if (@entries) {
    $populations{ $name } = $entries[0];
    return;
  }
  else {
    die "Could not add $name to populations\n";
  }
}
