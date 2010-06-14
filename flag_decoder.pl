#!/usr/bin/perl 
# 
# flag decoder code for SAM format
# 
# 
# Kim Brugger (07 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

while (<>) {
  chomp;

  print "-"x20 ."\n". value2text($_) . "-"x20 ."\n";
}


# 
# 
# 
# Kim Brugger (07 Apr 2010)
sub value2text {
  my ($value) = @_;

  my $text = "";

  if ($value & 0x0001) {
    $text .= "the read is paired in sequencing\n";
  }

  if ($value & 0x0002 ) {
    $text .= "the read is mapped in a proper pair\n";
  }
  if ($value & 0x0004 ) {
    $text .= "the query sequence itself is unmapped\n";
  }
  if ($value & 0x0008 ) {
    $text .= "the mate is unmapped\n";
  }
  if ($value & 0x0010 ) {
    $text .= "strand of the query (1 for reverse)\n";
  }
  if ($value & 0x0020 ) {
    $text .= "strand of the mate\n";
  }
  if ($value & 0x0040 ) {
    $text .= "the read is the first read in a pair\n";
  }
  if ($value & 0x0080 ) {
    $text .= "the read is the second read in a pair\n";
  }
  if ($value & 0x0100 ) {
    $text .= "the alignment is not primary\n";
  }
  if ($value & 0x0200 ) {
    $text .= "the read fails platform/vendor quality checks\n";
  }
  if ($value & 0x0400 ) {
    $text .= "the read is either a PCR or an optical duplicate\n";
  }


  return $text;
}



__END__

