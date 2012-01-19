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
#use warnings;
use Getopt::Std;


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
use EASIH::Mail;
use EASIH::DONE;

my $debug = 0;
#$debug = 1;

my %opts;
getopts('hvd:', \%opts);

if($opts{h})
{
    $0 =~ s/.*\///;
    print STDERR "\n\nDescription: Checks run complete (not the file Run.completed) status for illumina runs.\n";
    print STDERR "\nDetails: script for picking recent file in the relevant directory - looking for \"Copying logs to network run folder\" in Events.log file - alert people about run complete status.\n";
    print STDERR "\nUsage: $0 -h<help>\n\n\n";
    exit;
}

my $verbose = $opts{v} || 0;


my $to = 'bics@easih.ac.uk,lab@easih.ac.uk'; #global
$to = 'sri.deevi@easih.ac.uk,kim.brugger@easih.ac.uk'; #global
#$to = 'kim.brugger@easih.ac.uk'; #global
#my $dir = "/seqs/illumina3/";
#my $dir = "/seqs/babraham/";
my $RunDir; #global
my $rid;

my @dirs  = ('/seqs/illumina2/', 
	     '/seqs/illumina3/', 
	     '/seqs/illumina4/',
    );

@dirs = split(/,/, $opts{ 'd' }) if ( $opts{ 'd' } );
my $easih_toolbox = '/software/installed/easih-toolbox/';
$easih_toolbox = '/home/kb468/easih-toolbox/' if ( $debug );

foreach my $dir ( @dirs ) {
  opendir(DIR, "$dir");
  my @files = grep(!/^\.|\.log$/, sort readdir(DIR));
  closedir(DIR);
 
  while( $RunDir = shift @files) {
    
    print "RID:: $RunDir --> $rid\n" if ( $verbose );


#    next if ($RunDir !~ /111123/);
    
    my $runfolder   = ${dir}.$RunDir;
    my $eventfile   = "$runfolder/Events.log"; 
    my $RTAcomp     = "$runfolder/RTAComplete.txt";
    my $Intensities = "$runfolder/Data/Intensities"; 
    my $Basecalls   = "$runfolder/Data/Intensities/BaseCalls"; 

    my %Execute = ('01-BCL2QSEQ'  => {'01-command'    => "cd $Intensities; /software/bin/setupBclToQseq.py --no-eamss --in-place --overwrite -o BaseCalls -b BaseCalls", 
				  '02-prestatus'  => "BCL2QSEQ_SETUP", 
				  '03-failstatus' => "BCL2QSEQ_SETUP_FAILED", 
				  '04-poststatus' => "BCL2QSEQ_DONE"},
	       
	       '02-MAKE'      => {'01-command'    => "cd $Basecalls; make -j 8;", 
				  '02-prestatus'  => "MAKE_STARTED", 
				  '03-failstatus' => "MAKE_FAILED", 
				  '04-poststatus' => "MAKE_DONE"},
	       
	       '03-BCL2FQ'    => {'01-command'    => "$easih_toolbox/scripts/BaseCalls2fq.pl -di $Basecalls", 
				  '02-prestatus'  => "BCL2FQ_STARTED", 
				  '03-failstatus' => "BCL2FQ_FAILED", 
				  '04-poststatus' => "BCL2FQ_DONE"},
	       
	       '04-QC_Report' => {'01-command'    => "$easih_toolbox/scripts/QC_report.pl -p illumina -f ", # exclude $fqfile param until you start processing it
				  '02-prestatus'  => "QC_REPORT_STARTED", 
				  '03-failstatus' => "QC_REPORT_FAILED", 
				  '04-poststatus' => "QC_REPORT_DONE"},
		   
		   
	       '05-QC_DB' => {'01-command'    => "$easih_toolbox/DONE/QC_report.pl -p ILLUMINA -r ",
			      '02-prestatus'  => "QC_DB_STARTED", 
			      '03-failstatus' => "QC_DB_FAILED", 
			      '04-poststatus' => "QC_DB_DONE"},
    );

        
    if( -e $eventfile || -e $RTAcomp || -e "$runfolder/RTAComplete.txt") {

      my $checkstring = 1 if (-e $RTAcomp || -e "$runfolder/RTAComplete.txt");
      $checkstring = `grep -c "Copying logs to network run folder" $eventfile` 
	  if (! $checkstring);
      chomp( $checkstring );
      
      ### Update the status for old and unfinished runs ###
      if(!$checkstring) {
	print "$RunDir: Run not finished\n" if ($verbose);
	if ( @files ) {
	  $rid = EASIH::DONE::add_run($RunDir, 'ILLUMINA');
	  my $last_status = EASIH::DONE::fetch_latest_offloading_status($rid); 
	  print "$RunDir setting status to RUN_ABORTED\n" if ($verbose);
	  
	  RunStatus("RUN_ABORTED") if (! $last_status);
	  next;
	}
      }

      $rid = EASIH::DONE::add_run($RunDir, 'ILLUMINA');
      my $last_status = EASIH::DONE::fetch_latest_offloading_status($rid); 
      print "$RunDir last status: $last_status\n" if ($verbose && $last_status);
      
      next if ($last_status && $last_status ne "RETRY_OFFLOAD");

	my $failed_status = EASIH::DONE::fetch_offloading_failure( $rid ) if ( $last_status );
	my $fixed_failed_program = 0;

	### Start post processing for new runs ###
	if($checkstring)
	{
	  ### Run Complete ###
	  if ( ! $failed_status ) {
	    RunStatus("RUN_COMPLETED");
	    SendEmail("RUN_COMPLETED","",1);
	  }

            ### Run Programs: BCL2QSEQ, MAKE and BCL2FQ ###
	    
	    for my $program(sort keys %Execute)
	    {
		
		next, if($program =~ /QC_/);

		if ( $last_status && $last_status eq "RETRY_OFFLOAD" && ! $fixed_failed_program ) {
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

	  print "$failed_status $fixed_failed_program\n";
	    

	  if (! $last_status || 
	      $last_status ne "RETRY_OFFLOAD" ||
	      ($last_status eq "RETRY_OFFLOAD" && ! $fixed_failed_program && 
	       $failed_status eq $Execute{'04-QC_Report'}{'03-failstatus'})) {



	    ### QC_Report ###
	    my @fqfiles = EASIH::DONE::fetch_files_from_rid( $rid );

	    foreach my $fqfile(@fqfiles)
	    {
		
		my @GARstring;
		
		push @GARstring, $Execute{'04-QC_Report'}{'01-command'} . "$fqfile";
		push @GARstring, $Execute{'04-QC_Report'}{'02-prestatus'};
		push @GARstring, $Execute{'04-QC_Report'}{'03-failstatus'};
		push @GARstring, $Execute{'04-QC_Report'}{'04-poststatus'};
		
		next, if(!GrabAndRun(@GARstring));
	    }

	  }
	  
	  if (! $last_status || 
	      $last_status ne "RETRY_OFFLOAD" ||
	      ($last_status eq "RETRY_OFFLOAD" && ! $fixed_failed_program && 
	       $failed_status eq $Execute{'05-QC_DB'}{'03-failstatus'})) {


	    ### QC_Report ###
	    my %fqfiles = EASIH::DONE::fetch_files_from_rid_by_sample( $rid );
	    use Data::Dumper;
	    if (1) {
	      foreach my $sample (keys %fqfiles)  {
		
		my $param = " -1 $fqfiles{ $sample }[0] ";
		$param   .= " -2 $fqfiles{ $sample }[1] " if ($fqfiles{ $sample }[1]);
		
		my @GARstring;
		
		push @GARstring, $Execute{'05-QC_DB'}{'01-command'} . "$param";
		push @GARstring, $Execute{'05-QC_DB'}{'02-prestatus'};
		push @GARstring, $Execute{'05-QC_DB'}{'03-failstatus'};
		push @GARstring, $Execute{'05-QC_DB'}{'04-poststatus'};
		
		next if( !GrabAndRun(@GARstring) );
	      }
	    }
	  }
	  
	  
	  ### QC done email ###
	  $last_status = EASIH::DONE::fetch_latest_offloading_status( $rid );
	    
	  if($last_status eq "QC_DB_DONE") {
	    
	    ### Processing done status and email ###
	    my $end = "PROCESSING_DONE";
	    RunStatus($end);
	    SendEmail($end,$end,1);
	  }
	}
    }
  NEXT_RUNFOLDER:
  }
}

###########################################################
sub RunStatus
{
    my($status) = @_;
    
    EASIH::DONE::add_offloading_status($rid, 
				       "$status");
}


###########################################################
sub SendEmail
{
    my($subject, $message,$success) = @_;
 
    $subject = "ERROR:  $subject" if ( ! $success );

    $message = "Failed Command: $message", if(! $success);

    if($success && $message eq "PROCESSING_DONE")
    {
	#chomp(my $log = `~kb468/easih-toolbox/logistics/logistics_log.pl $CurrentRunDir`), 

	my $log = EASIH::DONE::run_log( $rid );
	$message = "$message \n\n $log";
    }

    my $body = "\n\n\t\t *****This is an automated email***** \n\n $message \n\n\t\t\t *****The End***** \n\n";
         
    
    EASIH::Mail::send($to, 
		      "[easih-done] $subject - $RunDir", 
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
	SendEmail($fail_status, $command,0);
	return 0;
    }
    
    RunStatus($post_status);
    
    return 1;
}
