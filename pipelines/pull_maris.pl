#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (13 Sep 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my $hpc_path = "/scratch/easih/MARIS/";
my $opts = 'p:P:';
my %opts;
getopts($opts, \%opts);

my $project_path = "$hpc_path/$opts{p}" if ( $opts{p} );
$project_path = "$opts{P}" if ( $opts{P} );

die "Please provide either a projects (-p) of a full HPC path (-P)\n" if ( ! $project_path );

my $tmp_dir = create_tmp_dir();

system "scp login.hpc.cam.ac.uk:$project_path/*.done $tmp_dir/";

print "scp \"login.hpc.cam.ac.uk:$project_path/*.done\" $tmp_dir/ \n";

print "$tmp_dir\n";

opendir(my $done_dir, "$tmp_dir");
while ( my $done_file =  readdir($done_dir)) {
  next if ( $done_file =~ /^\."/);

  open (my $done_fh, "$tmp_dir/$done_file") || die "Could not open '$tmp_dir/$done_file': $!\n";
  while (<$done_fh>) {
    chomp;
    system "mkdir /data/$opts{p}/MARIS/" if ( ! -e "/data/$opts{p}/MARIS/" );
    system "scp login.hpc.cam.ac.uk:$project_path/$_ /data/$opts{p}/MARIS/";
  }
  close ($done_fh);
}


# 
# 
# 
# Kim Brugger (13 Sep 2012)
sub create_tmp_dir {
  
  use File::Temp;
  my ($tmp_fh, $tmp_dir) = File::Temp::tempfile(DIR => "/tmp" );
  close ($tmp_fh);
  system "rm -f $tmp_dir";
  system "mkdir $tmp_dir";

  
  return $tmp_dir;
}

