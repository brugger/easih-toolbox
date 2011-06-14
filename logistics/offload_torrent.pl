#!/usr/bin/perl 
# 
# mgion01
# sudo apg-get install  libdbd-pg-perl libdbi-perl libnet-daemon-perl libplrpc-perl
#
# mgpc17
# sudo apt-get install postgresql-client-8.4 postgresql-client-common libdbd-pg-perl
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk


use strict;
use warnings;
use Data::Dumper;

use DBI;


# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment. Needs to be prior to the use of EASIH* modules.
BEGIN {
  my $path = $0;
  $path =~ s/^\.\///;
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
use EASIH::Logistics;
use EASIH::Mail;

#EASIH::Logistics::connect();

my $res_folder = "/results/analysis/";


my $dbi = DBI->connect("DBI:Pg:dbname=iondb;host=mgion01.medschl.cam.ac.uk", 'ion') || die "Could not connect to database: $DBI::errstr";
my $q = 'select "fastqLink", "sffLink" ';
$q .= 'from rundb_results results join rundb_experiment experiments on results.experiment_id = experiments.id where "status" = \'Completed\'';
my $sth = $dbi->prepare( $q );
$sth->execute();

while (my @results = $sth->fetchrow_array()) {
  my ($fq_file, $sff_file ) = @results;
  
  # extracting the information we need for moving off the data.
  my $res_dir = $fq_file;
  $res_dir =~ s/(.*\/).*/$res_folder$1/;
  my $run_folder = $res_dir;
  $run_folder =~ s/.*\/(.*)\//$1/;
  
  $fq_file =~ s/.*\/(.*)/$1/;
  $sff_file =~ s/.*\/(.*)/$1/;

  next if ($fq_file !~ /Z03/);

  my $org_fq = $fq_file;

  # misnaming problem... Just for code development
  $fq_file  =~ s/Z03[-_]/Z03/g;


  my $last_status = EASIH::Logistics::fetch_latest_run_folder_status("$run_folder");
  # a cron script sorts out the permissions on a hourly basis, so if we failed there before lets try again.
  if (  $last_status && $last_status ne "WRONG_PERMISSIONS" ) {
    print STDERR "Runfolder: $run_folder last know status is: $last_status, skipping it\n";
    next;
  }

  # so far things looks good. Update the tracking data base with the new run being finished
  EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "FINISHED_RUN") if (!$last_status);

  # extract the sample name from the fq file
  my ($sample_name) = $fq_file =~ /_([A-Z]\d{6,7})_/;

  # no name? send an email, and skip to the next runfolder.
  if ( ! $sample_name ) {
    EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "WRONGLY_NAMED");
    
    EASIH::Mail::send('kim.brugger@easih.ac.uk', 
		      "[easih-data] Failed Ion Torrent offload ($run_folder)", 
		      "Cannot find the EASIH-sample name in run folder $run_folder with the output name: $fq_file\n");
    
    next;
  }

  my ($outfile, $error) = EASIH::Logistics::sample2outfilename("$sample_name");

  if ($error) {
    EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "WRONG_PERMISSIONS");
    
    EASIH::Mail::send('kim.brugger@easih.ac.uk', 
		      "[easih-data] Failed Ion Torrent offload ($run_folder)", 
		      "$error \n");
    next;
  }

  $outfile .=".fq";


  $fq_file = $org_fq;

  EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "START_OFFLOADING");
  if ( system "scp ionadmin\@mgion01:$res_dir/$fq_file $outfile 2> /dev/null" ) {
    EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "OFFLOADING_FAILED");
    EASIH::Mail::send('kim.brugger@easih.ac.uk', 
		      "[easih-data] Offloading error Ion Torrent run ($run_folder)", 
		      "failed command: 'scp ionadmin\@mgion01:$res_dir/$fq_file $outfile'\n");
    next;
  }
  EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "FINISHED_OFFLOADING");

  EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "START_COMPRESSION");
  if ( system "gzip $outfile") {
    EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "COMPRESSION_FAILED");
    EASIH::Mail::send('kim.brugger@easih.ac.uk', 
		      "[easih-data] Data compression error on Ion Torrent data ($run_folder)", 
		      "failed command: 'gzip $outfile'\n");
    next;
  }
  EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "FINISHED_COMPRESSION");

  EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "START_QC");
  if ( system "/software/installed/easih-toolbox/scripts/QC_report.pl -p TORRENT -rf $outfile.gz") {
    EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "QC_FAILED");
    EASIH::Mail::send('kim.brugger@easih.ac.uk', 
		      "[easih-data] QC report failed on Ion Torrent data ($run_folder)", 
		      "failed command: '/software/installed/easih-toolbox/scripts/QC_report.pl -p TORRENT -rf $outfile.gz'\n");
    next;
  }
  EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "FINISHED_QC");

  EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "OFFLOADED");

  EASIH::Mail::send('kim.brugger@easih.ac.uk', 
		    "[easih-data] Successfully offloaded Ion Torrent data ($run_folder)", 
		    "$fq_file offloaded to $outfile.gz\n");

  EASIH::Logistics::add_file_offload($run_folder, $fq_file, "$outfile.gz");
}
