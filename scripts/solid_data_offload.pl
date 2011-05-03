#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (03 May 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts('cl:h', \%opts);
my $live_run      = $opts{ n } || 0;
my $compress_file = $opts{ c } || 0;

my $in_dir = shift;

if ( $opts{h}) {
  $0 =~ s/.*\///;
  print STDERR "\n";
  print STDERR "USAGE: script for copying data from a solid result folder to the shared disc. The script \n";
  print STDERR "USAGE: validates file names, and ensure that a counter is added if a this is a rerun\n";
  print STDERR "USAGE: $0 -c<ompress> -l<ive run, actually copies> -h<help> [input dir, default = latest result directory]\n";
  print STDERR "\n";
  exit -1;
}

my $hostname = `hostname`;
chomp $hostname;

if ( $hostname !~ /solid\d+\.sequencer/) {
  die "$hostname is not a SOLiD headnode\n";
}

$hostname =~ s/\.sequencer//;


if ( ! $in_dir ) {
  
   $in_dir ="/data/results/$hostname/". `ls -rt  /data/results/$hostname/ | tail -1`;
  chomp $in_dir;
}


open (my $find, "find $in_dir | ") ||
    die "Could not open find pipe: $!\n";


my @files;

while (<$find>) {
  chomp;
  next if ( !/csfasta\z|qual\z/);
  next if (/missing|unassigned/);
  next if (!/F3|F5|R3|R5/);

  my $org_file = $_;
  my ($path, $file) = ("", $_);
  ($path, $file) = ($1, $2)  if ( /(.*)\/(.*)/);

  print "$path/$file ";

  if ( $file !~ /.*?([A-Z]\d{7}_\d+).*/ && $file !~ /.*?([A-Z]\d{7}).*/) {
    print STDERR "Wrongly named file: $file\n";
    next;
  }


  #solid0354_20110421_FRAG_BC_A16_2_A16_2_F3_A1600009.csfasta
  if ( $file =~ /([R|F][3|5])_QV/) {
    my $read_type = $1;
    $file =~ s/.*?([A-Z]\d{7}_\d+).*/$1/ || $file =~ s/.*?([A-Z]\d{7}).*/$1/;
    $file .= "_$read_type\_QV.qual";
  }
  elsif ( $file =~ /([R|F][3|5]).csfasta/ || 
	  $file =~ /([R|F][3|5])_[A-Z]\d{7}.csfasta/ ||
	  $file =~ /([R|F][3|5])_[A-Z]\d{7}_\d+.csfasta/
      ) {
    my $read_type = $1;
    $file =~ s/.*?([A-Z]\d{7}_\d+).*/$1/ || $file =~ s/.*?([A-Z]\d{7}).*/$1/;
    $file .= "_$read_type\.csfasta";
  }

  my ($project) = $file =~ /([A-Z]\d{2})/;
  
 RENAMED_FILE:


  if ( -e "/ifs/data/$project/raw/$file" || -e "/ifs/data/$project/raw/$file.gz" ) {
    my ($sample, $version, $postfix) = $file =~ /([A-Z]\d{7})_(\d+)_(.*)/;
    ($sample, $version, $postfix) = $file =~  /([A-Z]\d{7})_(.*)/;
    
    if (! $postfix) {
      $postfix = $version;
      $version = 2;
    }

    $file = join("_", $sample, $version, $postfix);
    goto RENAMED_FILE;
  }  

  print "--> /ifs/data/$project/raw/$file\n" if ( ! $compress_file );
  print "--> /ifs/data/$project/raw/$file.gz\n" if ( $compress_file );

  if ( $live_run ) {
    system "mkdir -p /ifs/data/$project/raw/";
    system "cp $org_file /ifs/data/$project/raw/$file";
    system "gzip /ifs/data/$project/raw/$file" if ( $compress_file );
  }
}
  

