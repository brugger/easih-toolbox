#!/usr/bin/perl 
# 
# Map/realign/indels/snps integrated pipeline
# 
# 
# Kim Brugger (27 Jul 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use Getopt::Std;

use lib '/home/kb468/easih-pipeline/modules';
use EASIH::JMS;
use EASIH::JMS::Misc;
use EASIH::JMS::Samtools;
use EASIH::JMS::Picard;



our %analysis = ('csfasta2fastq'    => { function   => 'csfasta2fastq',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=05:00:00"},

		 'fastq-split'      => { function   => 'fastq_split',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=02:00:00"},
		 
		 'std-aln'          => { function   => 'bwa_aln',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=12:00:00"},
		 
		 'std-generate'     => { function   => 'bwa_generate',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},

		 'std-tag_sam'          => { function   => 'sam_add_tags',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},
		 
		 'std-sam2bam'          => { function   => 'EASIH::JMS::Samtools::sam2bam',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=1000mb,walltime=10:00:00"},
		 
		 'std-merge'     => { function   => 'EASIH::JMS::Picard::merge',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00",
					 sync       => 1},

		 'get_mapped'   => { function   => 'EASIH::JMS::Samtools::get_mapped',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'get_unmapped' => { function   => 'EASIH::JMS::Samtools::get_unmapped',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 

		 're-align'      => { function   => 'bwa_aln',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=12:00:00"},
		 
		 're-generate'     => { function   => 'bwa_generate',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},

		 're-sam2bam'          => { function   => 'EASIH::JMS::Samtools::sam2bam',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=1000mb,walltime=10:00:00"},
		 
		 'initial-merge'     => { function   => 'EASIH::JMS::Picard::merge',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00",
					 sync       => 1},

		 'initial-sort'         => { function   => 'EASIH::JMS::Picard::sort',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=20000mb,walltime=08:00:00"},
		 

		 'initial_index'    =>{ function   => 'EASIH::JMS::Samtools::index',
					hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
		 



		 'get_all_mapped'   => { function   => 'EASIH::JMS::Samtools::get_mapped',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'get_all_unmapped' => { function   => 'EASIH::JMS::Samtools::get_unmapped',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'seq_names'        => { function   => 'fetch_seq_names'},
		 
		 'identify_indel'   => { function   => 'identify_indel',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'realign_indel'    => { function   => 'realign_indel',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},

		 'merge_realigned'     => { function   => 'EASIH::JMS::Picard::merge',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00",
					 sync       => 1},

		 'sort_realigned'    => { function   => 'EASIH::JMS::Picard::sort',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=20000mb,walltime=08:00:00"},
		 

		 'index_realigned'   =>{ function   => 'EASIH::JMS::Samtools::index',
					hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
		 
		 'call_indels'      => { function   => 'call_indels',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},

		 'merge_indels'     => { function   => 'merge_indels',
					 sync       => 1},


		 'identify_snps'    => { function   => 'identify_snps',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'filter_snps'      => { function   => 'filter_snps',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500b,walltime=01:00:00"},
		 

		 'merge_vcfs'       => { function   => 'merge_vcfs',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=02:00:00", 
					 sync       => 1},
		 
		 'cluster_snps'     => { function   => 'cluster_snps',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=01:00:00"},
		 
		 'rescore_snps'     => { function   => 'rescore_snps',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=50000mb,walltime=01:00:00"},		 
		 
		 

    );
		     


our %flow = ( 'csfasta2fastq'    => "std-aln",
	      'fastq-split'      => "std-aln",

	      'std-aln'          => "std-generate",
	      'std-generate'     => "std-tag_sam",
	      'std-tag_sam'      => 'std-sam2bam',
	      'std-sam2bam'      => "std-merge",
	      'std-merge'        => ["get-mapped", "get-unmapped"],

	      "get-mapped"   => 'initial_merge',
	      "get-unmapped" => 're-aln',
	      're-aln'           => 're-generate',
	      're-generate'      => 're-sam2bam',
	      're-sam2bam'       => 'initial_merge',

	      'initial_sort'         => "initial_index",

	      "initial_index"    => ['get_all_mapped', 'get_all_unmapped'],
	      'get_all_unmapped' => 'merge_realigned',
	      'get_all_mapped'   => 'seq_names',
	      'seq_names'        => 'identify_indel',
	      'identify_indel'   => 'realign_indel',
	      'realign_indel'    => 'merge_realigned',
	      'merge_realigned'  => 'sort_realigned',
	      'sort_realigned'   => 'index_realigned',

	      'index_realigned'  => ['call_indels','identify_snps'],
	      'call_indels'      => 'merge_indels',

	      'identify_snps'    => 'filter_snps',
	      'filter_snps'      => 'merge_vcfs',
	      'merge_vcfs'       => 'cluster_snps',
	      'cluster_snps'     => 'rescore_snps'
	      );

EASIH::JMS::print_flow('fastq-split');

my %opts;
getopts('i:b:f:n:hlr:g:a:m:a:', \%opts);

if ( $opts{ r} ) {
  EASIH::JMS::restore_state($opts{r});
  getopts('i:b:f:n:hlr:g:', \%opts);
}


my $infile      = $opts{'i'} || $opts{'g'} || usage();
my $outfile     = $opts{'b'} || usage();
my $reference   = $opts{'f'} || usage();
my $split       = $opts{'n'} || 10000000;
my $align_param = $opts{'a'} || " ";

my $readgroup = $opts{'r'};
my $platform  = uc($opts{'p'}) || usage();


$platform = "SOLEXA" if ( $platform eq "ILLUMINA");
$align_param .= " -c " if ( $platform eq "SOLID");


my $bwa          = EASIH::JMS::Misc::find_program('bwa');
my $fq_split     = EASIH::JMS::Misc::find_program('fastq_split.pl');
my $samtools     = EASIH::JMS::Misc::find_program('samtools');
my $solid2fq     = EASIH::JMS::Misc::find_program('solid2fastq.pl');
my $tag_sam      = EASIH::JMS::Misc::find_program('tag_sam.pl');

#my $solid2fq  = '/home/kb468/bin/solid2fastq.pl';
#my $bwa       = '/home/easih/bin/bwa';
#my $samtools  = '/home/easih/bin/samtools';
#my $fq_split  = '/home/kb468/bin/fastq_split.pl';


#EASIH::JMS::verbosity(10);

#EASIH::JMS::dry_run('csfasta2fastq');
#exit;
EASIH::JMS::hive('Darwin');

if ( $opts{'g'} ) {
  
  &EASIH::JMS::run('BWA-mapping');
#  &EASIH::JMS::run('SAM2BAM');
  
}
else {
  &EASIH::JMS::run('csfasta2fastq');
}

&EASIH::JMS::store_state();

my $extra_report = "infile ==> $infile\n";
$extra_report .= "outfile ==> $outfile\n";
$extra_report .= "align_param ==> $align_param\n";
$extra_report .= "Binaries used..\n";
$extra_report .= `ls -l $samtools`;
$extra_report .= `ls -l $bwa` . "\n";

EASIH::JMS::mail_report('kim.brugger@easih.ac.uk', $outfile, $extra_report);



sub csfasta2fastq {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file();
 
  my $cmd = "$solid2fq";
  $cmd .= " -n $split " if ( $split );
  $cmd .= "$infile $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}


sub fastq_split {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file();
 
  my $cmd = "$solid2fq";
  $cmd .= " -n $split " if ( $split );
  $cmd .= "$infile $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}



sub bwa_aln {
  my ($input) = @_;
  
  $input = $opts{'g'} if ($opts{'g'});
  $input =~ s/\.\z//;

  my @inputs = glob "$input*";

  return if ( ! @inputs );

  foreach my $input ( @inputs ) {
    my $tmp_file = EASIH::JMS::tmp_file(".sai");
    my $cmd = "$bwa aln $align_param  -f $tmp_file $reference $input ";
    EASIH::JMS::submit_job($cmd, "$tmp_file $input");
  }
}

sub bwa_generate {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".sam");
  my $cmd = "$bwa samse -f $tmp_file $reference $input  ";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

# 
# Adds readgroup & aligner tags to the sam-file.
# 
# Kim Brugger (22 Jul 2010)
sub sam_add_tags {
  my ($input) = @_;

  if ( ! $readgroup ) {
    $readgroup = $infile;
    $readgroup =~ s/.fastq//;
    $readgroup =~ s/.fq//;
    $readgroup =~ s/.gz//;
  }

  my $cmd = "$tag_sam -R $input -p platform -r $readgroup -a bwa -A '$align_param' ";
  EASIH::JMS::submit_job($cmd, $input);
}


sub rename {
  my ($input) = @_;

  EASIH::JMS::submit_system_job("mv $input $outfile", $outfile);
}

# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {
  
  print "Not the right usage, please look at the code\n";
  exit;

}
