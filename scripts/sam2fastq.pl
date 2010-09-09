#!/usr/bin/perl
# Produce a fastq file from a SAM file

use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts('2:o:h', \%opts);

my $VERSION = '20091221.01';

sub process;

process;

exit;

# ----------------------------------------------------------------------
sub process {

  my @fhw = ();

  
  if ( $opts{o} || $opts{2}) {
    
    if ( ! $opts{2} ) {
      open( $fhw[0], "> $opts{o}") || die "Could not write to '$opts{o}': $!\n";
    }
    else {
      $opts{o} ||= $opts{2};
      open( $fhw[0], "> $opts{o}.1") || die "Could not write to '$opts{o}': $!\n";
      open( $fhw[1], "> $opts{o}.2") || die "Could not write to '$opts{o}': $!\n";
    }
  }
      

  while (my $line = <>) {    # sets $line

    chomp $line;

    next if $line =~ /^\@/xms;           # skip header lines

    my ($name, $flag, $read, $qual) = (split /\t/, $line)[0,1,9,10];

    if ($flag & 0x0010) {                # SAM revcomps reads aligned to reverse strand
      $read = reverse $read;
      $read =~ tr/ACGTacgt/TGCAtgca/;
      $qual = reverse $qual;
    }

    if ($flag & 0x0080 ) {
      $name .="/2";
    }
    else {
      $name .="/1";
    }

    if ( $opts{2} ) {
      # second read
      if ($flag & 0x0080 ) {
	print {$fhw[1]} "\@$name\n$read\n+\n$qual\n" ;
      }
      else {
	print {$fhw[0]} "\@$name\n$read\n+\n$qual\n" ;
      }

    }
    elsif ( $opts{o}) {
      print {$fhw[0]} "\@$name\n$read\n+\n$qual\n";
    }
    else {
      print  "\@$name\n$read\n+\n$qual\n";
    }


  }

  if ( $opts{2} ) {
    print STDOUT "$opts{o}.1\t$opts{o}.2\n" if (! -z "$opts{o}.2");
    print STDOUT "$opts{o}.1\n"             if ( -z "$opts{o}.2");
  }
  elsif ($opts{o}) {
    print STDOUT  "$opts{o}\n";
  }

  return;

}
