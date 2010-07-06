#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (07 Jun 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my %known;
my %unknown;

my @infos;

my $printed_header = 0;

while (<>) {


  chomp;

  push @infos, $1 if ( /##INFO=(.*?),/);

  next if (/\#/);

  my @fields = split("\t", $_);
  
  my ($chr, $pos, $id, $score, $stats ) = ($fields[0], $fields[1], $fields[2], $fields[5], $fields[7]);

  next if ( ! $stats );



  if ( ! $printed_header ) {
    @infos = grep (!/AB|MQ0/, @infos);
    print join("\t","Pos", "dbsnp", "QUAL", @infos) . "\n";
    $printed_header++;
  }

  my $known_snp = 0;
  $id = "unknown" if ($id eq ".");
  $id = "known"   if ($id ne "unknown");

  my %line;

  print "$chr:$pos\t$id\t$score\t";

  foreach my $value ( split(";", $stats)) {
    my ( $key, $value) = split("=", $value);
    $line{ $key} =  $value;
  }  
  

  foreach my $field ( @infos ) {

    $line{ $field } = "0" if ( $field eq "HRun" && ! $line{$field} );

    $line{ $field } = "NA" if ( ! defined $line{$field} );
    print "$line{ $field }\t";
  }
  print "\n";

}
