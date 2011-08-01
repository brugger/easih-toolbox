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
use EASIH::DONE;
use EASIH::Sample;


my $infile = shift || die "Needs an infile\n";

my @path = split("/", $infile);
my $run_folder = $path[3];
my $run_id = EASIH::DONE::add_run($run_folder, 'ILLUMINA');

my %fid_cache = ();

print "Rundir: $run_folder --> $run_id\n";

my %multiplex_counter;

my $read2 = 0;

open (my $in, "$infile") || die "Could not open '$infile': $!\n";
while( <$in> ) {
  chomp;
  my @f = split(/\t/ );
  @f = split(/\s{3,}/ ) if ( !@f || @f == 1);

  my ($lane, $read_nr) = $f[0] =~ /lane (\d).(\d)/;
  $read2 = 1 if ( $read_nr == 2 );
  my ($sample, $total_reads, $pass_filter) = ($f[1],$f[2],$f[3]);
  if ( ! $read_nr) {

    my $fid = find_or_create_fid("$sample.1.fq");
    my ($lane) = $f[0] =~ /lane (\d)/;
    print "FID == $fid $sample.1.fq $lane\n";
    my (undef, $sample, $bcode, $count, $perc) = @f;
    $multiplex_counter{$sample} += $count;
    $perc =~ s/[\% ]//;
    
    EASIH::DONE::add_illumina_multiplex_stats( $run_id, $fid, $lane, 1, $sample, $bcode, $total_reads, undef, $perc);
    EASIH::DONE::add_illumina_multiplex_stats( $run_id, $fid, $lane, 2, $sample, $bcode, $total_reads, undef, $perc) if ( $read2 );

  }
  else {
    print "$lane, $read_nr\n";
    $pass_filter =~ s/^(\d+).*/$1/;
    
    my $project = $sample =~ /^(\w{3})/;
    print "$sample.1.fq, $sample, $project, $run_folder, 'ILLUMINA', $total_reads, $pass_filter\n";
    
    my $fid = find_or_create_fid("$sample.$read_nr.fq");
    
    print "$run_id, $fid, $lane, $read_nr\n";
    
    EASIH::DONE::add_illumina_lane_stats( $run_id, $fid, $lane, $read_nr, $sample, $total_reads, $pass_filter )
	
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

  $filename = "/data/$project/raw/$filename.gz";


  my $fid = EASIH::DONE::add_file($filename, $sample, $project, $run_folder, 'ILLUMINA');
  $fid_cache{ $filename } = $fid;
  return $fid;
}

