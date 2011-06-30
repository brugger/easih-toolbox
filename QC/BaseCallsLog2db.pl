#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (24 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;



# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment. Needs to be prior to the use of EASIH* modules.
BEGIN {
  my $path = $0;
  if ($path =~ /.*\//) {
    $path =~ s/(.*)\/.*/$1/;
    push @INC, "$path/modules" if ( -e "$path/modules");
    $path =~ s/(.*)\/.*/$1/;
    push @INC, "$path/modules" if ( -e "$path/modules");
    
  }
  else {
    push @INC, "../modules" if ( -e "../modules");
    push @INC, "./modules" if ( -e "./modules");
  }
}


use EASIH;
use EASIH::QC::db;
use EASIH::Sample;


my $infile = shift || die "Needs an infile\n";

my @path = split("/", $infile);
my $run_folder = $path[3];
my $run_id = EASIH::QC::db::add_run($run_folder, 'ILLUMINA');

my %fid_cache = ();

print "Rundir: $run_folder --> $run_id\n";

my %multiplex_counter;

open (my $in, "$infile") || die "Could not open '$infile': $!\n";
while( <$in> ) {
  chomp;
  my @f = split("\t");
  
  my ($lane, $read_nr) = $f[0] =~ /lane (\d).(\d)/;
  my ($sample, $total_reads, $pass_filter) = ($f[1],$f[2],$f[3]);
  if ( ! $read_nr) {

    my $fid = find_or_create_fid("$sample.1.fq");
    my ($lane) = $f[0] =~ /lane (\d)/;
    print "FID == $fid $sample.1.fq $lane\n";
    my (undef, $sample, $bcode, $count, $perc) = @f;
    $multiplex_counter{$sample} += $count;
    $perc =~ s/[\% ]//;
    
    EASIH::QC::db::add_illumina_multiplex_stats( $run_id, $fid, $lane, $sample, $bcode, $total_reads, undef, $perc);

  }
  else {
    print "$lane, $read_nr\n";
    $pass_filter =~ s/^(\d+).*/$1/;
    
    my $project = $sample =~ /^(\w{3})/;
    print "$sample.1.fq, $sample, $project, $run_folder, 'ILLUMINA', $total_reads, $pass_filter\n";
    
    my $fid = find_or_create_fid("$sample.1.fq");
    
    print "$run_id, $fid, $lane, $read_nr\n";
    
    EASIH::QC::db::add_illumina_lane_stats( $run_id, $fid, $lane, $read_nr, $sample, $total_reads, $pass_filter )
	
  }

}




# 
# 
# 
# Kim Brugger (27 Jun 2011)
sub find_or_create_fid {
  my ($filename ) = @_;

  $filename =~ s/.*\///;
  print "$filename\n";

  return $fid_cache{ $filename } if ( $fid_cache{ $filename });

  my ($sample, $project) = EASIH::Sample::filename2sampleNproject($filename);
  
  if (! $sample ) {
    return undef;
  }


  my $fid = EASIH::QC::db::add_file($filename, $sample, $project, $run_folder, 'ILLUMINA');
  $fid_cache{ $filename } = $fid;
  return $fid;
}

