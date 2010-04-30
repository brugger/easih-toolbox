#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use Getopt::Std;

use lib '/home/kb468/projects/easih-flow/modules';
use lib '/home/kb468/easih-flow/modules';
use EASIH::JMS;
use EASIH::JMS::Samtools;

our %analysis = ('solid2bfast_fastq'   => { function   => 'solid2bfast_fastq',
					    hpc_param  => "-NEP-fqs -l mem=500mb,walltime=10:00:00"},
		 
		 'bfast-mapping'       => { function   => 'bfast_mapping',
					    hpc_param  => "-NEP-fqs -l mem=2000mb,walltime=03:00:00"},
		 
		 'bfast-localalign'    => { function   => 'bfast_localalign',
					    hpc_param  => "-NEP-fqs -l mem=2000mb,walltime=03:00:00"},
		 
		 'bfast-post'          => { function   => 'bfast_post',
					    hpc_param  => "-NEP-fqs -l mem=2000mb,walltime=03:00:00"},
		 
		 'SAM2BAM'             => { function   => 'EASIH::JMS::Samtools::sam2bam',
					    hpc_param  => "-NEP-fqs -l mem=2000mb,walltime=04:00:00"},
		 
		 'samtools-merge'      => { function   => 'EASIH::JMS::Samtools::bwa_merge',
					    hpc_param  => "-NEP-fqs -l mem=2000mb,walltime=05:00:00"},
		 
		 'samtools-sort'       => { function   => 'EASIH::JMS::Samtools::samtools_sort',
					    hpc_param  => "-NEP-fqs -l mem=2500mb,walltime=08:00:00"},
		 
		 'bam-rename'          => { function   => 'rename'},
		 
		 'samtools-index'      => { function   => 'samtools_index',
					    hpc_param  => "-NEP-fqs -l mem=2000mb,walltime=08:00:00"},

		 'find_input'          => { function   => 'find_input'},
		 );
		     
our %flow = ( 'solid2bfast_fastq' => "bfast-mapping",
	      'find_input'        => "bfast-mapping",
	      'bfast-mapping'     => "bfast-localalign",
	      'bfast-localalign'  => "bfast-post",
	      'bfast-post'        => "SAM2BAM",
	      'SAM2BAM'           => "samtools-merge",
	      'samtools-merge'    => "samtools-sort",
	      'samtools-sort'     => "bam-rename",
	      'bam-rename'        => "samtools-index");

#EASIH::JMS::job_stats('1007155');

my %opts;
getopts('w:i:q:o:p:n:r', \%opts);

if ( $opts{ r} ) {
  &EASIH::JMS::restore_state('bfast_pipeline.pl.freeze');
#  EASIH::JMS::print_HPC_usage();
#  exit;
  getopts('w:i:q:o:p:n:r', \%opts);
  EASIH::JMS::reset();
#  EASIH::JMS::dry_run('fastq-split');
#  exit;
}

my $infile  = $opts{'i'} || usage();
my $qual    = $opts{'q'} || usage();
my $outfile = $opts{'o'} || usage();
my $prefix  = $opts{'p'} || usage();
my $split   = $opts{'n'} || 10000000;

my $samtools    = '/home/kb468/bin/samtools';
my $solid2fastq = '/home/kb468/bin/solid2fastq';
my $bfast       = '/home/kb468/bin/bfast';

#EASIH::JMS::verbosity(10);
#EASIH::JMS::validate_flow('solid2bfast_fastq');
#EASIH::JMS::dry_run('find_input');
#exit;
#EASIH::JMS::run_flow('solid2bfast_fastq');
EASIH::JMS::run_flow('find_input');
&EASIH::JMS::store_state();

EASIH::JMS::delete_tmp_files();
EASIH::JMS::delete_hpc_logs();



# 
# 
# 
# Kim Brugger (28 Apr 2010)
sub find_input {
  
  my @fastq_files = glob "tmp/8RnXuXdHs5.*";
  EASIH::JMS::push_input( @fastq_files );

}




sub solid2bfast_fastq {
  my ($hpc_param) = @_;

  my $tmp_file   = EASIH::JMS::tmp_file();
  my $out_prefix = EASIH::JMS::tmp_file("", 1);
  
  my $cmd = "$solid2fastq -n $split -Z -o $out_prefix $infile $qual > $tmp_file 2> /dev/null";

  EASIH::JMS::submit_n_wait_job($cmd, $hpc_param);


  my @fastq_files = glob "$out_prefix.*";


  EASIH::JMS::tag_for_deletion(@fastq_files);
  EASIH::JMS::push_input( @fastq_files );

}

# 
# 
# 
# Kim Brugger (27 Apr 2010)
sub bfast_mapping {
  my ($hpc_param) = @_;

  my @inputs = EASIH::JMS::fetch_n_reset_inputs();

  my @cmds;
  foreach my $input ( @inputs ) {
    push @cmds, "$bfast match -f $prefix -A1 -k 8 -M 100 -Q2500  -z -r $input > $input.bmf";
    EASIH::JMS::push_input("$input.bmf");
  }


  EASIH::JMS::submit_n_wait_jobs(\@cmds, $hpc_param);
}

# 
# 
# 
# Kim Brugger (27 Apr 2010)
sub bfast_localalign {
  my ($hpc_param) = @_;

  my @inputs = EASIH::JMS::fetch_n_reset_inputs();

  my @cmds;
  foreach my $input ( @inputs ) {
    push @cmds, "$bfast localalign -f $prefix -A1 -o20 -m $input > $input.baf";
    EASIH::JMS::push_input("$input.baf");
  }

  EASIH::JMS::submit_n_wait_jobs(\@cmds, $hpc_param);
}


# 
# 
# 
# Kim Brugger (27 Apr 2010)
sub bfast_post {
  my ($hpc_param) = @_;

  my @inputs = EASIH::JMS::fetch_n_reset_inputs();

  my @cmds;
  foreach my $input ( @inputs ) {
    push @cmds, "$bfast postprocess -f $prefix-Q1000  -z -i $input > $input.sam";
    EASIH::JMS::push_input("$input.sam");
  }

  EASIH::JMS::submit_n_wait_jobs(\@cmds, $hpc_param);
}

sub rename {
  my ($hpc_param) = @_;
  my @inputs = EASIH::JMS::fetch_n_reset_inputs();
  my $cmd = "mv @inputs $outfile ";
  eval { system "$cmd" };
  EASIH::JMS::fail($@) if ($@);
}

sub samtools_index {
  my ($hpc_param) = @_;
  my $cmd = "$samtools index $outfile";
  EASIH::JMS::submit_job($cmd, $hpc_param);
  EASIH::JMS::wait_jobs( );  
}		     



# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {
  
  print "Not the right usage, please look at the code\n";
  exit;

}
