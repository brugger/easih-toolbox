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
use EASIH::JMS::Misc;
use EASIH::JMS::Samtools;



our %analysis = ('seq_names'      => { function   => 'fetch_seq_names'},
		 
		 'identify_indel' => { function   => 'identify_indel',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},
		 
		 'realign_indel'  => { function   => 'realign_indel',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},
		 
		 'bam-merge'      => { function   => 'EASIH::JMS::Samtools::merge',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00", 
				       sync       => 1},
		 
		 'samtools-sort'  => { function   => 'EASIH::JMS::Samtools::sort',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=50000mb,walltime=08:00:00"},
		 
		 'bam-rename'     => { function   => 'rename'},
		 
		 'samtools-index' => { function   => 'EASIH::JMS::Samtools::index',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},

		 'call_indels'    => { function   => 'call_indels',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},
    );
		     

our %flow = ( 
	      'identify_indel' => "realign_indel",
	      "realign_indel"  => 'bam-merge',
	      'bam-merge'      => "samtools-sort",
	      'samtools-sort'  => "bam-rename",
	      'bam-rename'     => "samtools-index",
              "samtools-index" => 'call_indels');

my %opts;
getopts('b:B:R:ho:r:', \%opts);

if ( $opts{ r} ) {
  EASIH::JMS::restore_state($opts{r});
    getopts('b:B:R:ho:r:', \%opts);
}


my $bam_file   = $opts{'b'} || usage();
my $new_bam    = $opts{'B'} || usage();
my $reference  = $opts{'R'} || usage();
my $report     = $opts{'o'} || usage();

my $samtools  = EASIH::JMS::Misc::find_program('samtools');
my $gatk      = EASIH::JMS::Misc::find_program('gatk ');

EASIH::JMS::hive('Darwin');
EASIH::JMS::max_retry(0);

&EASIH::JMS::run('identify_indel');

&EASIH::JMS::store_state();


sub identify_indel {

  my @names = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }

  foreach my $name ( @names ) {
    my $tmp_file = EASIH::JMS::tmp_file(".intervals");
    my $cmd = "$gatk -T RealignerTargetCreator -R $reference -o $tmp_file -I $bam_file -L $name";
    EASIH::JMS::submit_job($cmd, "$tmp_file -L $name");
  }
  
}

sub realign_indel {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd = "$gatk -T IndelRealigner -targetIntervals $input --output $tmp_file -R $reference -I $bam_file";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}


sub bam_merge { 
  my (@inputs) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".merged.bam");

  print "MERGE :: @inputs \n";

  if (@inputs == 1 ) {
    EASIH::JMS::submit_job("mv @inputs $tmp_file", $tmp_file, 1);
  }
  else {
    my $cmd = "$samtools merge $tmp_file @inputs ";
    print "$cmd \n";
    EASIH::JMS::submit_job($cmd, $tmp_file, 1);
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

  EASIH::JMS::submit_system_job("mv $input $new_bam", undef, 1);
}

sub samtools_index {
  my ($input) = @_;

  my $cmd = "$samtools index $new_bam";
  EASIH::JMS::submit_job($cmd);
}		     

sub call_indels {
  my ($input) = @_;

  my $cmd = "$gatk -T IndelGenotyperV2 -R $reference -O $report -I $new_bam";
  EASIH::JMS::submit_job($cmd);
}		     

# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {

#  return;
  
  print "Not the right usage, please look at the code\n";
  exit;

}

