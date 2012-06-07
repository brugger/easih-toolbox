#!/usr/bin/perl 
# 
# reads snp information from the local database
# 
# 
# Kim Brugger (26 Nov 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

require DBI;
my $dbi = DBI->connect('DBI:mysql:dbsnp132', 'easih_admin', 'easih') || die "Could not connect to database: DBI::errstr";

my $rs = shift || die "please call with an rs id\n";

#print "select * from snp rs = '$rs'\n";

my $sth = $dbi->prepare("select * from snp where rs = '$rs'");
$sth->execute();

while (my $ref = $sth->fetchrow_hashref ) {

  print "$$ref{rs} --> $$ref{chr}:$$ref{pos} $$ref{ref_base}>$$ref{alt_base}\n";

  my $sth2 = $dbi->prepare("select * from population where rs='$rs';");
  $sth2->execute();

  my ( $observations, $seen, $freq, $pop) = (0,0, 0);
  
  while (my $ref2 = $sth2->fetchrow_hashref ) {
    $pop++;
    my $count = int($$ref2{allele_freq}*$$ref2{sample_size});
    $observations += $$ref2{sample_size};
    $seen += $count;
    $freq += $$ref2{allele_freq};
    print "$$ref2{pop} --> $$ref2{sample_size}/$$ref2{allele_freq}/$count\n";
  }

  printf("TOT ==> $observations/%.2f/$seen\n", $freq/$pop);

  foreach my $flag (split(";", $$ref{flags})) {
    
    my ($id, $value) = split("=", $flag);


    next if ( grep /$id/, ("RV", "dbSNPBuildID", "VP", "WGT", "SLO"));


    my $sth3 = $dbi->prepare("select full from flags where id='$id';");
    $sth3->execute();
    my $ref3 = $sth3->fetchrow_hashref();
    
    if ( $value ) {
      print "$$ref3{full} == $value\n";
    }
    else {
      print "$$ref3{full}\n";
    }

  }  

}
