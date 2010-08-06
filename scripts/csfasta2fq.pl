#!/usr/bin/perl -w

# Author: lh3
# Note: Ideally, this script should be written in C. It is a bit slow at present.
# Also note that this script is different from the one contained in MAQ.

# Changes by kb468
# new flag -n[ax entries] will split the output into multiple gz files.
#
# Somewhat optimized the read1 function so the code runs twice as fast
#
# The input is now more flexible
# 
# kb468: Hacked beyond recognition. Actually an almost rewrite of the program...

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my %opts;
getopts('p:o:nh', \%opts);

usage() if ($opts{h});

my $prefix   = $opts{p} || usage();
my $out      = $opts{o} || usage();
# we compress by default, -n turns that off
my $compress = $opts{n} ? 0 : 1;

# strip away the (F3/R3).* postfixes so we are sure to have the correct basename
$prefix =~ s/F3.*//;
$prefix =~ s/R3.*//;

my (@fhr, @fhw);
my @fn_suff = ('F3.csfasta', 'F3_QV.qual', 'R3.csfasta', 'R3_QV.qual');
my $is_paired = (-f "$prefix$fn_suff[2]" || -f "$prefix$fn_suff[2].gz")? 1 : 0;
if ($is_paired) { # paired end
  for (0 .. 3) {
	my $fn = "$prefix$fn_suff[$_]";
	$fn = "gzip -dc $fn.gz |" if (!-f $fn && -f "$fn.gz");
	open($fhr[$_], $fn) || die("** Fail to open '$fn'.\n");
  }
  
  # for bwa to work with csfasta files, the R3 should be read1 and F3 read2. This is 
  # counter intuitive, but that is Chinese logic for you!
  if ( $compress ) {
    open($fhw[0], "|gzip >$out.2.fastq.gz")  || die; # this is NOT a typo
    open($fhw[1], "|gzip >$out.1.fastq.gz")  || die;
  }
  else {
    open($fhw[0], " >$out.2.fastq")  || die; # this is NOT a typo
    open($fhw[1], " >$out.1.fastq")  || die;
  }
  
  my (@df, @dr);
  @df = &read1(1); 
  @dr = &read1(2);
  while (@df && @dr) {

    print {$fhw[0]} $df[1]; 
    print {$fhw[1]} $dr[1];
      
    @df = &read1(1); 
    @dr = &read1(2);
  }

  if (@df) {
    print {$fhw[0]} $df[1];
    while (@df = &read1(1, $fhr[0], $fhr[1])) {
      print {$fhw[0]} $df[1];
    }
  }
  if (@dr) {
    print {$fhw[1]} $dr[1];
    while (@dr = &read1(2, $fhr[2], $fhr[3])) {
      print {$fhw[1]} $dr[1];
    }
  }
  close($fhr[$_]) for (0 .. $#fhr);
  close($fhw[$_]) for (0 .. $#fhw);
} else { # single end

  for (0 .. 1) {
    my $fn = "$prefix$fn_suff[$_]";
    $fn = "gzip -dc $fn.gz |" if (!-f $fn && -f "$fn.gz");
    open($fhr[$_], $fn) || die("** Fail to open '$fn'.\n");
  }
  if ( $compress ) {
    open($fhw[2], "| gzip >$out.1.fastq.gz") || die;
  }
  else {
    open($fhw[2], " >$out.1.fastq") || die;
  }

  my @df;
  while (@df = read1(1, $fhr[0], $fhr[1])) {
    print {$fhw[2]} $df[1];
  }

  close($fhr[$_]) for (0 .. $#fhr);
  close($fhw[2]);
}

#
# optimized ky kb468
#
sub read1 {
  my $i = shift(@_);
  my $j = ($i-1)<<1;
  my ($key, $seq);
  my ($fhs, $fhq) = ($fhr[$j], $fhr[$j|1]);
  while (<$fhs>) {
	my $t = <$fhq>;
	if (/^>(\d+)_(\d+)_(\d+)_[FR]3/) {
	  $key = sprintf("%.4d_%.4d_%.4d", $1, $2, $3); # this line could be improved on 64-bit machines
	  die(qq/** unmatched read name: '$_' != '$t'\n/) unless ($_ eq $t);
	  my $name = "$out:$1_$2_$3/$i";	  
	  $_ = substr(<$fhs>, 2);
	  tr/0123./ACGTN/;
	  my $s = $_;
	  $_ = <$fhq>;
	  s/^(\d+)\s*//;
	  s/-1/0/g;
	  my $qual;
	  map { $qual .= chr($_+33)} split(/\s+/);
	  $seq = qq/\@$name\n$s+\n$qual\n/;
	  last;
	}
  }
  return defined($seq)? ($key, $seq) : ();
}




# 
# 
# 
# Kim Brugger (06 Aug 2010)
sub usage {

  $0 =~ s/.*\///;
  
  print "USAGE: $0 -p[refix of the infiles] -o[outfile prefix] -n[o compress]\n";
  exit -1;

}