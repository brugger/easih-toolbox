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
my $dbi = DBI->connect('DBI:mysql:dbsnp132', 'easih_admin', 'easih') || die "Could not connect to database: DBI::errstr";

my $pop;
my $ref;
my $pop_warning = 0;
while (<>) {

  if (/##reference=(\w+)/) {
    $ref = $1;
    $ref = check_ref( $ref );
    next;
  }
  if (/##dbSNP_LOC_POP_ID=HapMap-(\w+)/i) {
    $pop = $1;
    next;
  }
  

  next if (/#/);
  die "could not find eitherref ($ref)\n" if ( !$ref);

  if ( !$pop && !$pop_warning ) {
    print STDERR  "could not find pop flag, just storing the basic snp information.\n";
    $pop_warning++;
  }
  chomp;
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
  my ($chr, $pos, $id, $ref_base, $alt_base, undef, undef, $info, $format) = split("\t");

  my $multi_genotypes = 0;
  $multi_genotypes++ if ( $format && $format =~ /GT2/);


  $info =~ s/NS=\d+;//;
  $info =~ s/VP=.*?;//;
  $info =~ s/AF=\d+\.\d+;//;


  my $sth = $dbi->prepare("select rs from snp where rs='$id'");
#  print "select rs from snp where rs='$id'\n";
  $sth->execute();
  my @entries =  $sth->fetchrow_array();
  if (  int (@entries) == 0) {
    $sth = $dbi->prepare("INSERT INTO snp (rs, chr, pos, flags, ref_id, ref_base, alt_base, multi_genotypes) VALUES ('$id', '$chr', '$pos', '$info', '$ref', '$ref_base', '$alt_base', '$multi_genotypes');");
    $sth->execute() || die"\n";
  }

  if ( $pop ) {
    my $allele_freq = $1 if ($info  =~ /AF=(\d+\.\d+)/);
    my $sample_size = $1 if ($info  =~ /NS=(\d+)/);
    $sth = $dbi->prepare("INSERT INTO population (pop, rs, sample_size, allele_freq) VALUES ('$pop', '$id', '$sample_size', '$allele_freq');");
    $sth->execute() || die "$pop', '$id'\n";
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
