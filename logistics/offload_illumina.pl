#!/usr/bin/perl -w
#*************************************************************************
#
#   File:       offload_illumina.pl
#   
#   Version:    V1.0
#   Date:       16.06.11
#   Function:   Offloading illumina data - postprocessing - running QC reports.
#   
#   Copyright:  (c) University of Cambridge / S.V.V. Deevi 2011
#   Author:     S.V.V. Deevi
#   Address:    Eastern Sequence Informatics Hub (EASIH),
#               The University of Cambridge,
#               Level 4, Laboratory Block Box 268,
#               Addenbrooke's Hospital,
#               Cambridge,
#               CB2 0QQ,
#               UK.
#   EMail:      svvd2@cam.ac.uk (or) sri.deevi@easih.ac.uk
#               
#
#   Description:
#   ============
#   script picking runfolders and going through them - if new, start post processing - alert people about run processing failures and update statuses on database.
#
#*************************************************************************
#
#   Usage:
#   ======
#   $0 -h<help>
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  16.06.11 Original.
#*************************************************************************

use lib "/home/kb468/easih-toolbox/modules/";
use strict;
use Getopt::Std;
use EASIH;
use EASIH::Mail;
use EASIH::Logistics;

my %opts;
getopts('h', \%opts);

if($opts{h})
{
  $0 =~ s/.*\///;
  print STDERR "\n\nDescription: Checks run complete (not the file Run.completed) status for illumina runs.\n";
  print STDERR "\nDetails: script for picking recent file in the relevant directory - looking for \"Copying logs to network run folder\" in Events.log file - alert people about run complete status.\n";
  print STDERR "\nUsage: $0 -h<help>\n\n\n";
  exit;
}

my $to = 'sri.deevi@easih.ac.uk,kim.brugger@easih.ac.uk';
#my $to = 'bics@easih.ac.uk,lab@easih.ac.uk';
#my $dir = "/seqs/illumina2/";
my $dir = "/seqs/babraham/";
my $CurrentRunDir;


opendir(DIR, "$dir");
my @files = readdir(DIR);
closedir(DIR);

 
foreach my $file (@files) 
{
    next, if($file =~ /^\.|\.log$/);
    
    my $runfolder   = ${dir}.$file;
    my $eventfile   = ${runfolder}."/Events.log"; 
    my $Data        = ${runfolder}."/Data";
    my $Intensities = ${Data}."/Intensities"; 
    my $Basecalls   = $Intensities."/BaseCalls"; 
    $CurrentRunDir  = $file;

    
    if(-e "$eventfile")
    {
	my $last_status = EASIH::Logistics::fetch_latest_run_folder_status("$file"); 
	
	next, if($last_status);
	
	chomp(my $checkstring = `grep -c "Copying logs to network run folder" $eventfile`);
	
	if($checkstring)
	{
	    ### Run Complete ###
	    my $body = "\n\n\t\t *****This is an automated email***** \n\n Run Completed for Illumina Run: $runfolder \n\n\t\t\t *****The End***** \n\n";
	    
	    EASIH::Mail::send($to, "[easih-data]Run Completed - $file", "$body");
	    RunStatus("RUN_COMPLETED");


            ### BCL2QSEQ ###
	    next, if(!GrabAndRun("cd $Intensities; /software/bin/setupBclToQseq.py --in-place --overwrite -o BaseCalls -b BaseCalls", 
				 "BCL2QSEQ_SETUP", 
				 "BCL2QSEQ_SETUP_FAILED", 
				 "BCL2QSEQ_DONE"));

	    
            ### MAKE ###
	    next, if(!GrabAndRun("cd $Basecalls; make;",
				 "MAKE_STARTED", 
				 "MAKE_FAILED",
				 "MAKE_DONE"));
	    

	    ### BCL2FQ ###
	    next, if(!GrabAndRun("/software/installed/easih-toolbox/scripts/BaseCalls2fq.pl -di $Basecalls > $Basecalls/BaseCalls2fq.log",
				 "BCL2QFQ_STARTED", 
				 "BCL2QFQ_FAILED",
				 "BCL2QFQ_DONE"));
	

	    ### Something Bad ###
	    my $last_status = EASIH::Logistics::fetch_latest_run_folder_status($file);
	    
	    if($last_status =~ /WRONG|FAIL/)
	    {
		next;
	    }
	    

	    
	    ### QC_Report ###
	    my @fqfiles = EASIH::Logistics::fetch_files_from_rundir($file);

	    foreach my $fqfile(@fqfiles)
	    {
		next, if(!GrabAndRun("/software/installed/easih-toolbox/scripts/QC_report.pl -p illumina -f $fqfile.1.fq.gz  -r",
				     "QC_REPORT_STARTED", 
				     "QC_REPORT_FAILED",
				     "QC_REPORT_DONE"));
	    }
	}
    }
}


###########################################################
sub RunStatus
{
    my($status) = @_;
    
    EASIH::Logistics::add_run_folder_status($CurrentRunDir, "ILLUMINA", "$status");
}


###########################################################
sub FailEmail
{
    my($command) = @_;
    
    EASIH::Mail::send($to, "[easih-data] Illumina Data Processing Failed!", "Failed Command: '$command'\n");
}


###########################################################
sub GrabAndRun
{
    my ($command, $pre_status, $fail_status, $post_status) = @_;
    
    RunStatus($pre_status);
    
    if(system($command))
    {
	RunStatus("$fail_status");
	FailEmail($command);
	return 0;
    }
    
    RunStatus($post_status);
    
    return 1;
}
