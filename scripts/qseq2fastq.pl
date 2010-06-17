#!/usr/bin/perl 

# Convert a qseq file to a fastq file. Assumes scores in qseq are
# phred+64. Optionally output only PF or non-PF reads.

use strict;
use warnings;
use English qw(-no_match_vars);
use Carp;
use Getopt::Long;

my $VERSION = '20100113.01';


sub process;
sub usage;
sub initialise;
sub timetag;

my $opts = initialise;
my $test = $opts->{test};

process;

exit;

# ----------------------------------------------------------------------
sub usage {

  ## no critic

  print STDERR "\n";
  print STDERR "swift_extract.pl version $VERSION\n";
  print STDERR "\n";
  print STDERR "    options:\n";
  print STDERR "\n";
  print STDERR "    --filter       filter out non-PF reads\n";
  print STDERR "    --no-filter    output only non-PF reads\n";
  print STDERR "    --run          use this run id in read names\n";
  print STDERR "    --verbose      include step-by-step stats\n";
  print STDERR "    --help         print this message and quit\n";
  print STDERR "\n";
  print STDERR "    Output to stdout, Default is to output all reads.\n";
  print STDERR "\n";

}

# ----------------------------------------------------------------------
sub process {

  my $log10 = log(10);
  my $count_in  = 0;
  my $count_out = 0;
  my %qtab;

  my $do_filter = exists $opts->{'filter'};
  my $no_filter = exists $opts->{'no-filter'};

  printf STDERR "%s%d qseq2fastq starting...\n", timetag, $PID, $count_in, $count_out;

  while (my $line = <>) {    # sets $line

    chomp $line;
    ++$count_in;

    my ($instr, $run_id, $lane, $tile, $x, $y, $index, $read, $bases, $q_line, $filter)
        = split /\t/, $line;

    if (exists $opts->{'run'}) {   # force run id? 
      $run_id = $opts->{'run'};
    }

    if ( ($filter and ! $no_filter) or ( ! $filter and ! $do_filter) ) {

      $bases =~ tr/./N/;           # turn dots into Ns

      ##### Commented out stuff below does logs-odds to phred conversion.
      ####$q_line =~ tr/!-\175/!!!!!!!!!!!!!!!!!!!!!!""""""##$$%%&&-++,-\136/;
      $q_line =~ tr/!-\175/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-\136/;

      if ($index ne '0') {
        print "\@${instr}_$run_id:$lane:$tile:$x:$y\#$index/$read\n";
      } else {
        print "\@${instr}_$run_id:$lane:$tile:$x:$y/$read\n";
      }
      print $bases, "\n";
      print "+\n";
      print $q_line, "\n";

      ++$count_out;

    }

  }

  printf STDERR "%s%d qseq2fastq read %d reads, wrote %d\n", timetag, $PID, $count_in, $count_out;

}

# ----------------------------------------------------------------------
sub initialise {

  my %opts;
  my $rc = GetOptions(\%opts, 'help', 'verbose', 'filter', 'no-filter', 'run=i');
  if ( ! $rc) {
    print {*STDERR} "\nerror in command line parameters\n" or croak 'print failed';
    usage;
    exit 1;
  }

  if (exists $opts{'help'}) {
    usage;
    exit;
  }

  if (exists $opts{'filter'} and exists $opts{'no-filter'}) {
    print {*STDERR} "\n--filter and --no-filter are mutually exclusive.\n\n";
    usage;
    exit 1;
  }

  return \%opts;

}

# ----------------------------------------------------------------------
sub timetag {

  my ($s, $m, $h, $D, $M, $Y) = localtime (time());
  return (sprintf "%04d%02d%02d-%02d%02d%02d ", $Y+1900, $M+1, $D, $h, $m, $s);

}
