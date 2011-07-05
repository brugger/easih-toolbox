#!/usr/bin/perl 
# 
# mgion01
# sudo apg-get install  libdbd-pg-perl libdbi-perl libnet-daemon-perl libplrpc-perl
#
# mgpc17
# sudo apt-get install postgresql-client-8.4 postgresql-client-common libdbd-pg-perl
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk
# Edited by Sri (16 Jun 2011), contact: sri.deevi@easih.ac.uk


use strict;
use warnings;
use Data::Dumper;

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
my $DYNAMIC_LIB_PATHS = 1;
BEGIN {
  if ( $DYNAMIC_LIB_PATHS ) {
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
  else {
    use lib '/home/kb468/easih-toolbox/modules/';
  }
}

use EASIH;
use EASIH::Logistics;
use EASIH::Sample;
use EASIH::Mail;

my $res_folder = "/results/analysis/";
my $CurrentRunDir; # Global ### svvd 
my $to = 'kim.brugger@easih.ac.uk'; #global
my $mailcount; #global

### Fetch Data from iondb (located on Ion torrent) ###
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
  
  $Currentrundir = $run_folder;
  $mailcount = 0;

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
  if(!$last_status)
  {
      my $body = "\n\n\t\t *****This is an automated email***** \n\n Run Completed for Ion torrent Run: $runfolder \n\n\t\t\t *****The End***** \n\n";
      SendEmail("$body"); 
      RunStatus("RUN_COMPLETED");
  } ### svvd 

  # extract the sample name from the fq file
  my ($sample_name) = $fq_file =~ /_([A-Z]\d{6,7})_/;
  
  # no name? send an email, and skip to the next runfolder.
  if ( ! $sample_name ) {
      RunStatus("WRONGLY_NAMED"); ### svvd
      
      SendEmail("Cannot find the EASIH-sample name in run folder $run_folder with the output name: $fq_file\n"); ### svvd
    
      next;
  }

  my ($outfile, $error) = EASIH::Sample::sample2outfilename("$sample_name");

  if ($error) {
    RunStatus("WRONG_PERMISSIONS"); ### svvd
    
    SendEmail("$error"); ### svvd
    next;
  }

  $outfile .=".fq";


  $fq_file = $org_fq;

  next, if(!GrabAndRun("scp ionadmin\@mgion01:$res_dir/$fq_file $outfile 2> /dev/null",
		       "OFFLOADING_STARTED",
		       "OFFLOADING_FAILED",
		       "OFFLOADING_DONE")); ### svvd
  
  
  next, if(!GrabAndRun("gzip $outfile",
		       "COMPRESSION_STARTED",
		       "COMPRESSION_FAILED",
		       "COMPRESSION_DONE")); ### svvd

  EASIH::Logistics::add_file_offload($run_folder, $fq_file, "$outfile.gz");
  my @fqfiles = EASIH::Logistics::fetch_files_from_rundir($file);

  foreach my $fqfile(@fqfiles)
  {
      next, if(!GrabAndRun("/software/installed/easih-toolbox/scripts/QC_report.pl -p illumina -f $fqfile -r",
			   "QC_REPORT_STARTED", 
			   "QC_REPORT_FAILED",
			   "QC_REPORT_DONE"));
  }

  #EASIH::Logistics::add_run_folder_status($run_folder, "TORRENT", "OFFLOADED");

  #EASIH::Mail::send('kim.brugger@easih.ac.uk', 
		    "[easih-data] Successfully offloaded Ion Torrent data ($run_folder)", 
		    "$fq_file offloaded to $outfile.gz\n");

}


### svvd ###
###########################################################
sub RunStatus
{
    my($status) = @_;

    EASIH::Logistics::add_run_folder_status($CurrentRunDir, 
					    "TORRENT", 
					    "$status");
}


### svvd ###
###########################################################
sub SendEmail
{
    my($message) = @_;
    
    if($mailcount)
    {
	EASIH::Mail::send($to, 
			  "[easih-data] Ion Torrent Data Processing Failed - $CurrentRunDir!", 
			  "Failed Command: '$message'\n");
    }
    else
    {
	EASIH::Mail::send($to, 
			  "[easih-data] Ion Torrent Run Completed - $CurrentRunDir", 
			  "$message");
    }
    
    $mailcount++;  
}


### svvd ###
###########################################################
sub GrabAndRun
{
    my ($command, $pre_status, $fail_status, $post_status) = @_;
    
    RunStatus($pre_status);
    
    if(system($command))
    {
	RunStatus("$fail_status");
	SendEmail($command);
	return 0;
    }
    
    RunStatus($post_status);
    
    return 1;
}
