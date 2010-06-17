#!/usr/bin/perl 
# 
# splits a fastq file into smaller sub files
# 
# 
# Kim Brugger (22 Apr 2010), contact: kim.brugger@easih.ac.uk


use strict;
use warnings;
use Getopt::Long;


my $entries = 10000000; # 10 mill.
my $prefix;

my $hm = GetOptions('entries:n'   => \$entries,
		    'prefix:s'    => \$prefix);

my $infile = shift || usage();
$prefix = $infile . "_split" if ( ! $prefix);
my $file_counter  = 1;
my $entry_counter = 0;


open (my $fastq, " $infile") || die "Could not open '$infile': $!\n";
open (my $outfile, "> $prefix." . $file_counter) || die "could not open '$prefix$file_counter': $!\n";
print "$prefix.$file_counter\n";
$file_counter++;
while (<$fastq>) {
  
  $_ .= <$fastq>;
  $_ .= <$fastq>;
  $_ .= <$fastq>;

  print $outfile $_;
  
  if ( ++$entry_counter >= $entries ) {
    close $outfile;
    open ( $outfile, "> $prefix." . $file_counter) || die "could not open '$prefix.$file_counter': $!\n";
    print "$prefix.$file_counter\n";
    $file_counter++;
    $entry_counter = 0;
  }
}
close( $fastq);
close $outfile;

sub usage {

  $0 =~ s/.*\///;
  print "USAGE: $0 --entries <number, def=10mill> --prefix <name, def = [infile]_split.XX] fastq-file\n";
  exit;
}
