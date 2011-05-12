#!/usr/bin/perl 
# 
# adds centimorgan information to selected snps.
# 
# 
# Kim Brugger (26 Nov 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

require DBI;
my $dbi = DBI->connect('DBI:mysql:dbsnp132', 'easih_admin', 'easih') || die "Could not connect to database: DBI::errstr";

my $pop;
my $ref = shift || die "Need a reference";
my $pop_warning = 0;
while (<>) {

  next if (/^\d+\t/);

  chomp;
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
  my ($chr, $rs,  undef, $cm) = split("\t");

  next if ( ! $cm);

  $chr = "X" if ( $chr == 23);

  my $sth = $dbi->prepare("select rs from snp where rs='$id'");
  $sth->execute();
  my @entries =  $sth->fetchrow_array();
  if (  int (@entries)) {
    
    $sth = $dbi->prepare("UPDATE snp set centimorgan='$cs';");
    $sth->execute() || die"\n";
  }

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
