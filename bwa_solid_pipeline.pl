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
use lib '/home/kb468/easih-pipeline/modules';
use EASIH::JMS;



our %analysis = ('csfasta2fastq'   => { function   => 'csfasta2fastq',
					hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=05:00:00"},
		 
		 'BWA-mapping'    => { function   => 'bwa_aln',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=03:00:00"},
		 
		 'BWA-samse'      => { function   => 'bwa_samse',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=25000mb,walltime=03:00:00"},
		 
		 'SAM2BAM'        => { function   => 'sam2bam',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=1000mb,walltime=03:00:00"},
		 
		 'BAM-merge'      => { function   => 'bam_merge',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},

#		 'remap'          => { function   => 'bam_remap',
#				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'samtools-sort'  => { function   => 'samtools_sort',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=50000mb,walltime=08:00:00"},
		 
		 'bam-rename'     => { function   => 'rename'},
		 
		 'samtools-index' =>{ function   => 'samtools_index',
				      hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
    );
		     

our %flow = ( 'csfasta2fastq' => "BWA-mapping",
	      'BWA-mapping'   => "BWA-samse",
	      'BWA-samse'     => "SAM2BAM",
	      'SAM2BAM'       => "BAM-merge",
#	      'BAM-merge'     => 'remap'
	      'BAM-merge'     => "samtools-sort",
	      'samtools-sort' => "bam-rename",
	      'bam-rename'    => "samtools-index");

my %opts;
getopts('i:b:f:n:hlr', \%opts);

if ( $opts{ r} ) {
  EASIH::JMS::restore_state('jms_test.pl.freeze');
#  EASIH::JMS::print_HPC_usage();
#  exit;
  getopts('w:i:o:p:n:r', \%opts);
#  EASIH::JMS::reset();
#  EASIH::JMS::dry_run('fastq-split');
#  exit;
}


my $infile    = $opts{'i'} || usage();
my $outfile   = $opts{'b'} || usage();
my $reference = $opts{'f'} || usage();
my $split     = $opts{'n'} || 30000000;

my $solid2fq  = '/home/kb468/bin/solid2fastq.pl';
my $bwa       = '/home/kb468/bin/bwa';
my $samtools  = '/home/kb468/bin/samtools';


#EASIH::JMS::verbosity(10);

#EASIH::JMS::dry_run('csfasta2fastq');
#exit;
EASIH::JMS::hive('Darwin');
EASIH::JMS::run('csfasta2fastq');
&EASIH::JMS::store_state();

EASIH::JMS::delete_tmp_files();
EASIH::JMS::delete_hpc_logs();


sub csfasta2fastq {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file();
  my $basename = $infile;
  
  $basename =~ s/F3.csfasta//;
  $basename =~ s/F3_QV.qual//;
  
  my $cmd = "$solid2fq";
  $cmd .= " -n $split " if ( $split );
  $cmd .= "$basename $tmp_file";

#  print STDERR "$cmd \n";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}



sub bwa_aln {
  my ($input) = @_;
  
  my @inputs = glob "$input.*";

  return if ( ! @inputs );

  foreach my $input ( @inputs ) {
    my $tmp_file = EASIH::JMS::tmp_file(".sai");
    my $cmd = "$bwa aln -e5 -c -f $tmp_file $reference $input ";
    EASIH::JMS::submit_job($cmd, "$tmp_file $input");
  }
}

sub bwa_samse {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".sam");
  my $cmd = "$bwa samse -f $tmp_file $reference $input  ";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub sam2bam {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd = "$samtools view -b -S $input > $tmp_file ";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub bam_merge { 
  my (@inputs) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".merged.bam");

  if (@inputs == 1 ) {
    system "mv @inputs $tmp_file";
  }
  else {
    my $cmd = "$samtools merge $tmp_file @inputs ";
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}

sub samtools_sort {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".merged.sorted");
  my $cmd = "$samtools sort -m 2048000000 $input $tmp_file ";
    
  EASIH::JMS::submit_job($cmd, "$tmp_file.bam");
}

sub rename {
  my ($input) = @_;

  my $cmd = "mv $input $outfile ";
  eval { system "$cmd" };
  EASIH::JMS::fail($@) if ($@);
}

sub samtools_index {
  my ($input) = @_;

  my $cmd = "$samtools index $outfile";
  EASIH::JMS::submit_job($cmd);
}		     

# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {
  
  print "Not the right usage, please look at the code\n";
  exit;

}
