#!/usr/bin/perl 
# 
# reads in the vcf header into the database
# 
# 
# Kim Brugger (26 Nov 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

require DBI;
my $dbi = DBI->connect('DBI:mysql:dbsnp132', 'easih_admin', 'easih') || die "Could not connect to database: DBI::errstr";


while (<>) {
  if ( /##INFO=\<ID=(\w+?),.*Description=\"(.*?)\"/ ) {
    my ($id, $full) = ($1, $2);

    my $sth = $dbi->prepare("INSERT INTO flags (id, full) VALUES (\"$id\",\"$full\");");
    $sth->execute();
  }
}

