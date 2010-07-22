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

# nohup ./scripts/indel_pipeline.pl -b 1.bam -B demo_01.realigned.bam -o 1.indels -R /home/easih/refs/hg18/hg18.fasta > realign.log &



our %analysis = ('seq_names'      => { function   => 'fetch_seq_names'},
		 
		 'dump_unmapped'  => { function   => 'dump_unmapped',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'identify_indel' => { function   => 'identify_indel',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'realign_indel'  => { function   => 'realign_indel',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'bam-merge'      => { function   => 'EASIH::JMS::Picard::merge',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=01:00:00", 
				       sync       => 1},
		 
		 'samtools-sort'  => { function   => 'EASIH::JMS::Picard::sort',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=50000mb,walltime=02:00:00"},
		 
		 'bam-rename'     => { function   => 'rename'},
		 
		 'samtools-index' => { function   => 'EASIH::JMS::Samtools::index',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},

		 'call_indels'    => { function   => 'call_indels',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},

		 'merge_indels'    => { function   => 'merge_indels',
					sync       => 1},
    );
		     

our %flow = ( 
	      "dump_unmapped"  => 'bam-merge',
	      'identify_indel' => "realign_indel",
	      "realign_indel"  => 'bam-merge',
	      'bam-merge'      => "samtools-sort",
	      'samtools-sort'  => "bam-rename",
	      'bam-rename'     => "samtools-index",
              "samtools-index" => 'call_indels',
	      'call_indels'    => 'merge_indels');



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

my $samtools    = EASIH::JMS::Misc::find_program('samtools');
my $gatk        = EASIH::JMS::Misc::find_program('gatk ');
my $fix_sorted  = EASIH::JMS::Misc::find_program('bam_fix_sorted_flag.pl');

EASIH::JMS::hive('Darwin');
EASIH::JMS::max_retry(0);

&EASIH::JMS::run('identify_indel', 'dump_unmapped');

&EASIH::JMS::store_state();

my $extra_report = "infile ==> $bam_file\n";
$extra_report .= "outfile ==> $new_bam\n";
$extra_report .= "report ==> $report\n";
$extra_report .= "Binaries used..\n";
$extra_report .= `ls -l $samtools`;
$extra_report .= `ls -l $gatk` . "\n";
$extra_report .= `ls -l $fix_sorted` . "\n";

EASIH::JMS::mail_report( 'kim.brugger@easih.ac.uk', "$report / $new_bam", $extra_report);

sub identify_indel {

  if ( ! -e $bam_file ) {
    die "$bam_file does not exist\n";
  }

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
    EASIH::JMS::submit_job($cmd, "$tmp_file $name");
  }
  
}

sub realign_indel {
  my ($input) = @_;

  my ($interval_file, $region) = split(" ", $input);

  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd;
  # If the interval file is empty the realigner ignores the region and produces an empty bamfile...
  if (  -z $interval_file ) {
    $cmd = "$samtools view -b $bam_file $region > $tmp_file";
  }
  else {
    $cmd = "$gatk -T IndelRealigner -targetIntervals $interval_file -L $region --output $tmp_file -R $reference -I $bam_file";
  }

  EASIH::JMS::submit_job($cmd, $tmp_file);
}


sub dump_unmapped {
  
  my $tmp_file = EASIH::JMS::tmp_file(".unmapped.bam");
  my $cmd = "$samtools view -bf4 $bam_file > $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub rename {
  my ($input) = @_;

  EASIH::JMS::submit_system_job("mv $input $new_bam", $new_bam);
}

sub call_indels {
  my ($input) = @_;

  my @names = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }

  foreach my $name ( @names ) {

    my $tmp_file = EASIH::JMS::tmp_file(".indels");
    my $cmd = "$gatk -T IndelGenotyperV2 -R $reference -O $tmp_file -I $new_bam -L $name ";
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}		     



# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub merge_indels {
  my (@inputs) = @_;


  print "MERGE :: @inputs \n";

  my $cmd = "cat  @inputs > $report";
  EASIH::JMS::submit_system_job("cat  @inputs > $report");
  
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

