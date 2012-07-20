#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (01 Dec 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;


# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 1;
  if ( $DYNAMIC_LIB_PATHS ) {
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules" && ! grep /^$path\/modules/, @INC);
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    push @INC, '/home/kb468/easih-toolbox/modules/';
  }
}


use EASIH;
use EASIH::Illumina::Sample_sheet;
use EASIH::Sample;

my $debug = 0;
#$debug = 1;

my %opts;
getopts("f:s:ho:", \%opts);

usage() if ($opts{h});
my $fastq_file   = $opts{f} || usage();
my $outdir      = $opts{'o'};
$outdir = "/tmp/BC2FQ/" if ($debug);
my $sample_sheet = "sample_sheet.csv" if (-e "sample_sheet.csv");
$sample_sheet    =  $opts{s} if ($opts{s} && -e $opts{s});
usage() if (! $sample_sheet );


my ($res, $errors ) = EASIH::Illumina::Sample_sheet::readin( $sample_sheet );
die $errors if ( $errors );

my %bcodes;
my %fhs;
foreach my $lane ( keys %$res ) {
  foreach my $bcode ( keys %{$$res{$lane}} ) {
    my $sample_name = $$res{$lane}{$bcode};
    my ($base_filename, $error) = EASIH::Sample::sample2outfilename( "$sample_name.1", $outdir);
    
    my $fh;
    open ($fh, "| gzip -c > $base_filename.fq.gz") || fail( "Could not open '$base_filename': $!\n", "BASECALL2FQ_PATH_ERROR");
    $fhs{ "$sample_name.1" }  = $fh;
    $bcodes{$bcode} = "$sample_name.1";

  }
}

open (my $fqs, $fastq_file) || die "Could not open '$fastq_file': $!\n";
while(<$fqs>) {
  my $header     = $_;
  my $sequence   = <$fqs>;
  my $strand     = <$fqs>;
  my $quality    = <$fqs>;

  my ($bcode) = ($sequence =~ /^(.{10})/);
  next if (! $bcode || ! $bcodes{$bcode});
  my $fout = $fhs{$bcodes{$bcode}};
  print $fout "$header$sequence$strand$quality";
  

}


# 
# 
# 
# Kim Brugger (01 Dec 2011)
sub usage {
  die "BLA BLA BLA\n";
}
