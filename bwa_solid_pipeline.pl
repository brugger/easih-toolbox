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

use lib '/home/kb468/easih-pipeline/modules';
use EASIH::JMS;
use EASIH::JMS::Misc;
use EASIH::JMS::Samtools;
use EASIH::JMS::Picard;



our %analysis = ('csfasta2fastq'  => { function   => 'csfasta2fastq',
					hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=05:00:00"},

		 'BWA-mapping'    => { function   => 'bwa_aln',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=12:00:00"},
		 
		 'BWA-samse'      => { function   => 'bwa_samse',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},

		 'tag_sam'        => { function   => 'sam_add_tags',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},
		 
		 'SAM2BAM'        => { function   => 'EASIH::JMS::Samtools::sam2bam',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=1000mb,walltime=10:00:00"},
		 
		 'bam-merge'      => { function   => 'EASIH::JMS::Picard::merge',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00",
				       sync       => 1},

		 'bam-sort'       => { function   => 'EASIH::JMS::Picard::sort',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=20000mb,walltime=08:00:00"},
		 
		 'bam-rename'     => { function   => 'rename'},
		 
		 'bam-index'      =>{ function   => 'EASIH::JMS::Samtools::index',
				      hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
    );
		     


our %flow = ( 'csfasta2fastq' => "BWA-mapping",
	      'BWA-mapping'   => "BWA-samse",
	      'BWA-samse'     => "tag_sam",
	      'tag_sam'       => 'SAM2BAM',
	      'SAM2BAM'       => "bam-merge",
	      'bam-merge'     => "bam-sort",
	      'bam-sort     ' => "bam-rename",
	      'bam-rename'    => "bam-index");

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

sub bwa_samse {
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
