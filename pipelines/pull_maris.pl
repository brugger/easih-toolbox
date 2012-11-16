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

my $project_dir = find_project_dir( $opts{p} );

opendir(my $done_dir, "$tmp_dir");
while ( my $done_file =  readdir($done_dir)) {
  next if ( $done_file =~ /^\."/);

  open (my $done_fh, "$tmp_dir/$done_file") || die "Could not open '$tmp_dir/$done_file': $!\n";
  while (<$done_fh>) {
    chomp;
    system "mkdir $project_dir/MARIS/" if ( ! -e "$project_dir/MARIS/" );
    system "scp login.hpc.cam.ac.uk:$project_path/$_ $project_dir/MARIS/";
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


# 
# 
# 
# Kim Brugger (14 Nov 2012)
sub find_project_dir {
  my $project_id = shift;

  return "/data/$project_id/" if ( -e "/data/$project_id/");
  
  my $letter  = substr($project_id, 0, 1);
  my $number  = substr($project_id, 1, 2);  

  my $project_dir = sprintf("/data/$letter/%02d_%02d/$project_id/", (int($number/10)*10), (int($number/10)*10+9));
  return $project_dir if ( -e $project_dir );

  die "Could not resolve project dir for $project_id\n";
  
}
