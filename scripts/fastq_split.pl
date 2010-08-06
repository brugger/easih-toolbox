#!/usr/bin/perl 
# 
# splits a fastq file into smaller sub files
# 
# 
# Kim Brugger (22 Apr 2010), contact: kim.brugger@easih.ac.uk


use strict;
use warnings;
use Getopt::Std;


my %opts;
getopts('1:2:e:p:h');

my $entries     = $opts{'e'} || 10000000; # 10 mill.
my $first_file  = $opts{'1'} || usage();
my $second_file = $opts{'2'} || usage();
my $prefix;
usage() if ( $opts{h});

my $hm = GetOptions('e:n'   => \$entries,
		    'prefix:s'    => \$prefix);




$prefix = $first_file . ".split" if ( ! $prefix);
my $file_counter  = 1;
my $entry_counter = 0;


open (my $first,  " $first_file") || die "Could not open '$first_file': $!\n";
open (my $second, " $second_file") || die "Could not open '$second_file': $!\n" if ( $second_file );
open (my $outfile, "> $prefix." . $file_counter) || die "could not open '$prefix$file_counter': $!\n";
print "$prefix.$file_counter\n";
$file_counter++;
while (<$first>) {
  
  $_ .= <$first>;
  $_ .= <$first>;
  $_ .= <$first>;

  print $outfile $_;
  
  if ( ++$entry_counter >= $entries ) {
    close $outfile;
    open ( $outfile, "> $prefix." . $file_counter) || die "could not open '$prefix.$file_counter': $!\n";
    print "$prefix.$file_counter\n";
    $file_counter++;
    $entry_counter = 0;
  }
}
close( $first);
close( $second) if ( $second_file);
close $outfile;

sub usage {

  $0 =~ s/.*\///;
  print "USAGE: $0 --entries <number, def=10mill> --prefix <name, def = [infile]_split.XX] -1 fq-file -2 fq-file\n";
  exit;
}
