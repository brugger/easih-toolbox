#!/usr/bin/perl 
# 
# Trims fq entries to the requested length.
# 
# 
# Kim Brugger (15 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH;

my %opts;
getopts('hiIvf:sro:', \%opts); 

Usage() if ($opts{"h"});

my $outfile  = $opts{o} || "";;
my $pattern = shift if (!$opts{f});
if (!$opts{f} && !$pattern) {
  print STDERR "No patteren or patteren file provided\n";
  Usage();
}

my @patterns;
# the patterns is provided in a file. 
if ($opts{f}) {
  open INFILE, "$opts{f}" or die "Could not open '$opts{f}': $!\n";
  while (<INFILE>) {
    chomp;
    push @patterns, $_;
  }
  close INFILE;
}

open (*STDOUT, ">$outfile") || die "Could not write to '$outfile': $!\n" if ( $outfile );


LOOP:
if (@ARGV) {
    $STDIN = shift;
    open (STDIN) or die "Could not open file $STDIN: $!";
}


while(my $header = <> ) {
  chomp( $header );
  my $sequence   = <>;
  chomp( $sequence );
  my $strand     = <>;
  chomp( $strand );
  my $quality    = <>;
  chomp( $quality );
  

  die "fastq file is either not a fq-file or contains extra lines.\n"  
      if ( $header !~ /^\@/ || $strand !~ /^[+-]/ );

  my $comp = $header;
  $comp = $sequence if ( $opts{s});
  
  my $comres = compare($comp);

  if ( !$comres && $opts{s} && $opts{r}) {
    $comp = revDNA($comp);
    $comres = compare($comp);
  }

  if (($opts{"v"})) {
    if (not $comres) {
      print join("\n", $header, $sequence, $strand, $quality). "\n";
    }
  }
  else  {
    if ( $comres) {
      print join("\n", $header, $sequence, $strand, $quality). "\n";
    }
  }

  

}   
close STDIN or die "could not close file";
goto LOOP if (@ARGV);



# 
# 
# 
# Kim Brugger (17 May 2011)
sub revDNA {
  my ($dna) = @_;

  $dna =~ tr/[ACGT]/[TGCA]/;
  $dna = reverse($dna);
  return $dna;
}



# 
# 
# 
# Kim Brugger (15 Jun 2011)
sub compare {
  my ($input) = @_;

  my $comres;
  # There is not one pattern but multiple ones
  if ($opts{f}) {
    foreach my $pattern (@patterns) {
    if ($opts{"i"}) {
      $comres = $input =~ /$pattern/i;
      }	
      elsif ($opts{"I"}) {
	$comres = $input =~ /^$pattern\z/i;
      }	
      else {
	$comres = $input =~ /$pattern/;
      }	
      last if $comres;
    }
  }
  else {
    if ($opts{"i"}) {
      $comres = $input =~ /$pattern/i;
    }	
    elsif ($opts{"I"}) {
      $comres = $input =~ /^$pattern\z/i;
    }	
    else {
      $comres = $input =~ /$pattern/;
    }	
  }

  return $comres;
}







# 
# 
# 
# Kim Brugger (16 May 2011)
sub Usage {

  $0 =~ s/.*\///;

  print "USAGE: $0 [patteren] infile(s) or stdin";
  print "USAGE: $0 greps either in the header or the sequence of a fastq file\n";
  print "USAGE: $0 -f<ile with patterens> \n";
  print "USAGE: $0 -h<elp, this message> \n";
  print "USAGE: $0 -i<gnore upper/lower case> \n";
  print "USAGE: $0 -o<ut file to write the results to> \n";
  print "USAGE: $0 -v< inverted match> \n";
  print "USAGE: $0 -I<D or sequence, perfect matches only> \n";
  print "USAGE: $0 -s<equence search> \n";
  print "USAGE: $0 -r<evert sequence for searching with both strands> \n";
  exit -1;
  
}

