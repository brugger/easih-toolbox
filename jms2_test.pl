#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use File::Temp;

use Getopt::Std;

my %job_stats;

my %opts;
getopts('w:i:o:p:n:', \%opts);

my $infile  = $opts{'i'} || "1_1_1.fastq" || usage();
my $outfile = $opts{'o'} || "1_1_1.fastq.bam" || usage();
my $prefix  = $opts{'p'} || " /home/kb468/Hs/hs_GRCh37" || usage();
my $split   = $opts{'n'} || 10000 || 30000000;

my $fq_split = '/home/kb468/bin/fastq_split.pl';
my $bwa      = '/home/kb468/bin/bwa';
my $samtools = '/home/kb468/bin/samtools';

use lib '/home/kb468/projects/easih-flow/modules';
use lib '/home/kb468/easih-flow/modules';
use EASIH::JMS;

our %analysis = ('fastq-split'   => { function   => 'fastq_split',
				      hpc_param  => "-NEP-fqs -l mem=500mb,walltime=00:05:00"},
		 
		 'BWA-mapping'    => { function   => 'bwa_aln',
				       hpc_param  => "-NEP-fqs -l mem=1000mb,walltime=00:10:00"},
		 
		 'BWA-samse'      => { function   => 'bwa_samse',
				       hpc_param  => "-NEP-fqs -l mem=1000mb,walltime=00:10:00"},
		 
		 'SAM2BAM'        => { function   => 'sam2bam',
				       hpc_param  => "-NEP-fqs -l mem=500mb,walltime=00:10:00"},
		 
		 'BWA-merge'      => { function   => 'bwa_merge',
				       hpc_param  => "-NEP-fqs -l mem=500mb,walltime=00:10:00",
				       sync       => 1},
		 
		 'samtools-sort'  => { function   => 'samtools_sort',
				       hpc_param  => "-NEP-fqs -l mem=1000mb,walltime=00:10:00"},
		 
		 'bam-rename'     => { function   => 'rename'},
		 
		 'samtools-index' =>{ function   => 'samtools_index',
				      hpc_param  => "-NEP-fqs -l mem=1000mb,walltime=00:10:00"},
    );
		     

our %flow = ( 'fastq-split'   => "BWA-mapping",
	      'BWA-mapping'   => "BWA-samse",
	      'BWA-samse'     => "SAM2BAM",
	      'SAM2BAM'       => "BWA-merge",
	      'BWA-merge'     => "samtools-sort",
	      'samtools-sort' => "bam-rename",
	      'bam-rename'    => "samtools-index");


#push @EASIH::JMS::jobs, '1003679';
#EASIH::JMS::wait_jobs( );

EASIH::JMS::verbosity(10);

#EASIH::JMS::validate_flow('fastq-split');
#exit;
EASIH::JMS::run('fastq-split');
#EASIH::JMS::dry_run('fastq-split');
#EASIH::JMS::delete_tmp_files();


sub fastq_split {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file();
  my $cmd = "$fq_split -e $split -p $tmp_file $infile";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub bwa_aln {
  my ($input) = @_;
  
  $input = " 1_1_1.fastq_split";
  
  my @inputs = glob "$input.*";

  return if ( ! @inputs );

  print "@inputs\n";

  foreach my $input ( @inputs ) {
    my $tmp_file = EASIH::JMS::tmp_file();
    my $cmd = "$bwa aln -f $tmp_file $prefix $input ";
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}

sub bwa_samse {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file();
  my $cmd = "$bwa samse $prefix $input > $tmp_file ";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub sam2bam {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd = "$samtools view -b -S $input > $tmp_file ";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub bwa_merge { 
  my (@inputs) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".merged.bam");
  my $cmd = "$samtools merge $tmp_file @inputs ";
  EASIH::JMS::submit_job($cmd, $tmp_file);
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
  
  print "not the right usage\n";
#  exit;

}
