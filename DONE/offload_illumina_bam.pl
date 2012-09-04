#!/usr/bin/perl -w
#*************************************************************************
#               
#
#   Description:
#   ============
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
$debug = 1;

my $easih_toolbox = '/software/installed/easih-toolbox/';
$easih_toolbox = '/home/kb468/easih-toolbox/' if ( $debug );


EASIH::DONE::Connect('done_dev') if ($debug); 

my %opts;
getopts('hvr:s:', \%opts);

usage() if ( $opts{h} );


use constant {
  STEP_NAME     => 0,
  PROGRAM       => 1
};

my @programs = (["BAM2BCL", "$easih_toolbox/scripts/bcl2bam.pl -dr "],
		["QC", "$easih_toolbox/DONE/QC_run.pl"]);


my $verbose = $opts{v} || 0;

my $to = 'bics@easih.ac.uk,lab@easih.ac.uk'; #global
$to = 'sri.deevi@easih.ac.uk,kim.brugger@easih.ac.uk'; #global
$to = 'kim.brugger@easih.ac.uk' if ( $debug );


my @sequencer_dirs  = ('/seqs/illumina2/', 
	     '/seqs/illumina3/', 
	     '/seqs/illumina4/',
	     '/seqs/illumina5/',
    );

# multiple specif dirs might have been given as input, rather than running as a deamon
@sequencer_dirs = split(/,/, $opts{ 's' }) if ( $opts{ 's' } );


foreach my $sequencer_dir ( @sequencer_dirs ) {

  opendir(my $sdir, "$sequencer_dir");
  # get an order list of the direcotories for this sequencer, ignoring the '.' and '..' dirs
  my @run_dirs = map {"$sequencer_dir/$_"} grep(!/^\./, sort readdir( $sdir ));
  closedir( $sdir );
  
  foreach my $run_dir ( @run_dirs ) {

    if (!is_an_illumina_run_dir( $run_dir )) {
#      print "$run_dir is not a full Illumin run directory \n";
      next;
    }

#    print "$run_dir is being investigated \n";

    if ( run_completed( $run_dir )) {
      print "$run_dir is completed\n";


      my $rid = EASIH::DONE::add_run($run_dir, 'ILLUMINA');
      my $last_status = EASIH::DONE::fetch_latest_offloading_status( $rid ); 
      $last_status ||= "Unknown";
      print "$run_dir last status: $last_status\n";


      next if ($last_status && $last_status ne "RETRY_OFFLOAD");

      failed_offloading() if (! run_program(@{$programs[0]}, $rid));
      failed_offloading() if (! run_program(@{$programs[1]}, $rid));
    }
  }
#  exit;
}


sub run_program {
  my ( $step_name, $command, $rid ) = @_;

  EASIH::DONE::add_offloading_status($rid, "$step_name\_START");

  
  print "$command $rid\n";
    
  if(system( "$command $rid"))   {
    EASIH::DONE::add_offloading_status($rid, "$step_name\_FAILED");
    send_failed_email($rid, "$step_name\_FAILED");
    return 0;
  }
    
  EASIH::DONE::add_offloading_status($rid, "$step_name\_DONE");
    
  return 1;
}



# 
# 
# 
# Kim Brugger (04 Sep 2012)
sub send_failed_email {
  SendEmail(@_, 0);
  
}

# 
# 
# 
# Kim Brugger (04 Sep 2012)
sub send_success_email {
  SendEmail(@_, 1);
  
}



###########################################################
sub SendEmail {
  return;
  my($rid, $subject, $message,$success) = @_;
 
  $subject = "ERROR:  $subject" if ( ! $success );
  
  $message = "Failed Command: $message", if(! $success);
  
  if($success && $message eq "PROCESSING_DONE") {
    #chomp(my $log = `~kb468/easih-toolbox/logistics/logistics_log.pl $CurrentRunDir`), 
    
    my $log = EASIH::DONE::run_log( $rid );
    $message = "$message \n\n $log";
  }
  
  my $body = "\n\n\t\t *****This is an automated email***** \n\n $message \n\n\t\t\t *****The End***** \n\n";
  
  my $run_dir = EASIH::DONE::fetch_run_name( $rid );

  
  EASIH::Mail::send($to, 
		    "[easih-done] $subject - $run_dir", 
		    "$body");
}


###########################################################


# 
# 
# 
# Kim Brugger (04 Sep 2012)
sub is_an_illumina_run_dir {
  my ( $run_folder ) = @_;

  $run_folder =~ s/\/{2,}/\//g;
  return 0 if ( ! -d $run_folder );
  return 0 if ( ! -e "$run_folder/Data/Intensities/BaseCalls");
  
  return 1;
}




# 
# 
# 
# Kim Brugger (04 Sep 2012)
sub run_completed {
  my ( $run_folder ) = @_;

  #This is for the HiSeq and MiSeq
  if ( -e "$run_folder/RTAComplete.txt" ) {
    return 1;
  }

  return 0 if ( ! -e "$run_folder/Events.log");
  
  # grep counting how many times the "Copying logs to network run folder" occurs on the log file: 1 == finished, 0 still running/aborted.
  my $checkstring = `grep -c "Copying logs to network run folder" $run_folder/Events.log`;
  chomp( $checkstring );
  return $checkstring;

}



# 
# 
# 
# Kim Brugger (04 Sep 2012)
sub usage {
  
  $0 =~ s/.*\///;
  print STDERR "\n\nDescription: Checks run complete (not the file Run.completed) status for illumina runs.\n";
  print STDERR "\nDetails: script for picking recent file in the relevant directory - looking for \"Copying logs to network run folder\" in Events.log file - alert people about run complete status.\n";
  print STDERR "\nUsage: $0 -h<help> -s<equencer directory> -D<irectory> -v<erbose>\n\n\n";
  exit -1;
}


__END__	   
