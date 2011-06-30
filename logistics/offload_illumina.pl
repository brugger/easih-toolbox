#!/usr/bin/perl -w
#*************************************************************************
#
#   File:       offload_illumina.pl
#   
#   Date:       22.06.11
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


use strict;
use Getopt::Std;


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


my $to = 'sri.deevi@easih.ac.uk,kim.brugger@easih.ac.uk'; #global
#my $to = 'bics@easih.ac.uk,lab@easih.ac.uk'; #global
my $dir = "/seqs/illumina2/";
#my $dir = "/seqs/babraham/";
my $CurrentRunDir; #global


opendir(DIR, "$dir");
my @files = grep(!/^\.|\.log$/, sort readdir(DIR));
closedir(DIR);

   
while(my $file = shift @files) 
{
    my $runfolder   = ${dir}.$file;
    my $eventfile   = ${runfolder}."/Events.log"; 
    my $Data        = ${runfolder}."/Data";
    my $Intensities = ${Data}."/Intensities"; 
    my $Basecalls   = $Intensities."/BaseCalls"; 
    $CurrentRunDir  = $file;    
    

    my %Execute = ('01-BCL2QSEQ'  => {'01-command'    => "cd $Intensities; /software/bin/setupBclToQseq.py --in-place --overwrite -o BaseCalls -b BaseCalls", 
				      '02-prestatus'  => "BCL2QSEQ_SETUP", 
				      '03-failstatus' => "BCL2QSEQ_SETUP_FAILED", 
				      '04-poststatus' => "BCL2QSEQ_DONE"},
		   
		   '02-MAKE'      => {'01-command'    => "cd $Basecalls; make;", 
				      '02-prestatus'  => "MAKE_STARTED", 
				      '03-failstatus' => "MAKE_FAILED", 
				      '04-poststatus' => "MAKE_DONE"},
		   
		   '03-BCL2FQ'    => {'01-command'    => "/home/kb468/easih-toolbox/scripts/BaseCalls2fq.pl -di $Basecalls > $Basecalls/BaseCalls2fq.log", 
				      '02-prestatus'  => "BCL2QFQ_STARTED", 
				      '03-failstatus' => "BCL2QFQ_FAILED", 
				      '04-poststatus' => "BCL2QFQ_DONE"},
		   
		   '04-QC_Report' => {'01-command'    => "/software/installed/easih-toolbox/scripts/QC_report.pl -p illumina -r -f ", # exclude $fqfile param until you start processing it
				      '02-prestatus'  => "QC_REPORT_STARTED", 
				      '03-failstatus' => "QC_REPORT_FAILED", 
				      '04-poststatus' => "QC_REPORT_DONE"},
	);

        
    if(-e "$eventfile")
    {
	my $last_status = EASIH::Logistics::fetch_latest_run_folder_status("$file"); 
	
	next, if($last_status && $last_status ne "RETRY_OFFLOAD");
	
	chomp(my $checkstring = `grep -c "Copying logs to network run folder" $eventfile`);
	

	### Update the status for old and unfinished runs ###
	if((!$checkstring) && (@files)) 
	{
	    RunStatus("RUN_ABORTED");
	    next;
	}


	my $failed_status = EASIH::Logistics::fetch_run_folder_failure( $file ) if ( $last_status );
	my $fixed_failed_program = 0;


	### Start post processing for new runs ###
	if($checkstring)
	{
	    ### Run Complete ###
	    RunStatus("RUN_COMPLETED") if ( ! $failed_status );
	    my $last_status = EASIH::Logistics::fetch_latest_run_folder_status($file);
	    SendEmail($last_status,1), if($last_status eq "RUN_COMPLETED");


            ### Run Programs: BCL2QSEQ, MAKE and BCL2FQ ###
	    
	    for my $program(sort keys %Execute)
	    {
		
		next, if($program =~ /QC_Report/);

		if ( $last_status eq "RETRY_OFFLOAD" && ! $fixed_failed_program ) {
		    next if ($failed_status ne $Execute{$program}{'03-failstatus'});
		    $fixed_failed_program++;
		}

		my @GARstring;
		
		for my $info(sort keys %{$Execute{$program}})
		{
		    push @GARstring, $Execute{$program}{$info};
		}
	
		goto NEXT_RUNFOLDER , if(!GrabAndRun(@GARstring));
	    }         
	    

	    ### Something Bad ###
	    $last_status = EASIH::Logistics::fetch_latest_run_folder_status($file);
	    #if($last_status =~ /WRONG|FAIL/)
	    if( $last_status ne "BASECALLS2FQ_DONE" )
	    {
		next;
	    }
	    

	    ### QC_Report ###
	    my @fqfiles = EASIH::Logistics::fetch_files_from_rundir($file);
	    
	    foreach my $fqfile(@fqfiles)
	    {
		
		my @GARstring;
		
		push @GARstring, $Execute{'04-QC_Report'}{'01-command'} . "$fqfile";
		push @GARstring, $Execute{'04-QC_Report'}{'02-prestatus'};
		push @GARstring, $Execute{'04-QC_Report'}{'03-failstatus'};
		push @GARstring, $Execute{'04-QC_Report'}{'04-poststatus'};
		
		next, if(!GrabAndRun(@GARstring));
	    }

    
	    
	    ### QC done email ###
	    $last_status = EASIH::Logistics::fetch_latest_run_folder_status($file);
	    
	    if($last_status eq "QC_REPORT_DONE")
	    {
		SendEmail($last_status,1), 
		
                ### Processing done status and email ###
		my $end = "PROCESSING_DONE";
		RunStatus($end);
		SendEmail($end,1);
	    }
	}
    }
  NEXT_RUNFOLDER:
}

###########################################################
sub RunStatus
{
    my($status) = @_;
    
    EASIH::Logistics::add_run_folder_status($CurrentRunDir, 
					    "ILLUMINA", 
					    "$status");
}


###########################################################
sub SendEmail
{
    my($message,$success) = @_;
 
    my $subject = $message;

    $message = "Failed Command: $message", if(! $success);

    if($success && $message eq "PROCESSING_DONE")
    {
	#chomp(my $log = `~kb468/easih-toolbox/logistics/logistics_log.pl $CurrentRunDir`), 

	my $log = EASIH::Logistics::runfolder_log($CurrentRunDir);
	$message = "$message \n\n $log";
    }

    my $body = "\n\n\t\t *****This is an automated email***** \n\n $message \n\n\t\t\t *****The End***** \n\n";
         
    
    EASIH::Mail::send($to, 
		      "[easih-data] Illumina Data Processing [Return Code = $success]: $subject - $CurrentRunDir", 
		      "$body");
}


###########################################################
sub GrabAndRun
{
    my ($command, $pre_status, $fail_status, $post_status) = @_;

    RunStatus($pre_status);
    
    if(system($command))
    {
	RunStatus("$fail_status");
	SendEmail($command,0);
	return 0;
    }
    
    RunStatus($post_status);
    
    return 1;
}
