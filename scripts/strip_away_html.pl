#!/usr/bin/perl 
# 
# Strips away html tags from report files, and transforms tables dividers into tabs
# 
# 
# Kim Brugger (17 Sep 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;


my %opts;
getopts('i:o:h', \%opts);
usage() if ( $opts{h});

my $infile    = $opts{i};
my $outfile   = $opts{o};

open (*STDIN, $infile) || die "Could not open '$infile': $!\n" if ( $infile );
open (*STDOUT, " > $outfile " ) || die "Could not open '$outfile': $!\n" if ( $outfile );

while (<>) {
  next if (/^<TABLE / || /^<\/TABLE /);
  
  $_ =~ s/<TR><TD>//;
  $_ =~ s/<\/T[RD]>//g;
  $_ =~ s/<TD>/\t/g;
  $_ =~ s/&nbsp;//g;
  $_ =~ s/<.*?>//g;
  
#  my @fields = split("<TR><TD>");

#  print join("\t", @fields);
  print;
}
