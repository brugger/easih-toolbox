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



our %analysis = ('call_indels'      => { function   => 'call_indels',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},

		 'merge_indels'     => { function   => 'merge_indels',
					 sync       => 1},
		 );
		     


our %flow = ( 'call_indels'      => 'merge_indels',
	      );

#EASIH::JMS::no_store();
#EASIH::JMS::print_flow('fastq-split');

my %opts;
getopts('b:R:o:p:h', \%opts);

# if ( $opts{ r} ) {
#   EASIH::JMS::restore_state($opts{r});
#   getopts('i:b:f:n:hlr:g:', \%opts);
# }


my $bam_file    = $opts{'b'} || usage();
my $reference   = $opts{'R'} || usage();
my $filters     = $opts{'f'} || "default";
my $report      = $opts{'o'} || usage();

my $readgroup   = $opts{'r'} || $report;
my $platform    = uc($opts{'p'}) || usage();
$platform = 'SOLEXA'      if ( $platform eq 'ILLUMINA');

# set platform specific bwa aln parameters
$align_param .= " -c "    if ( $platform eq "SOLID");
$align_param .= " -q 15 " if ( $platform eq "SOLEXA");

my $samtools     = EASIH::JMS::Misc::find_program('samtools');
my $gatk         = EASIH::JMS::Misc::find_program('gatk ');

validate_input();

#EASIH::JMS::verbosity(10);
EASIH::JMS::hive('Darwin');
#EASIH::JMS::hive('Kluster');
EASIH::JMS::max_retry(0);


&EASIH::JMS::run('call_indels');

&EASIH::JMS::store_state();

my $extra_report .= "bamfile ==> $bam_file\n";
$extra_report .= "indel_file ==> $report.indel\n";

$extra_report .= "align_param ==> $align_param + -e5 -t5 for second round aligning\n";
$extra_report .= "Binaries used..\n";
$extra_report .= `ls -l $samtools`;
$extra_report .= `ls -l $bwa` . "\n";

EASIH::JMS::mail_report('kim.brugger@easih.ac.uk', $bam_file, $extra_report);


sub fastq_split {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file();
 
  my $cmd = "$fq_split -e $split  -1 $first";
  $cmd .= " -2 $second "  if ( $second);
  $cmd .= " -o tmp/ > $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}



sub bwa_aln {
  my ($input) = @_;

  
  if ( $no_split ) {
    
    if ( $first && $second ) {
      my $first_tmp_file  = EASIH::JMS::tmp_file(".sai");
      my $second_tmp_file = EASIH::JMS::tmp_file(".sai");
      my $cmd = "$bwa aln $align_param  -f $first_tmp_file  $reference $first ;";
      $cmd   .= "$bwa aln $align_param  -f $second_tmp_file $reference $second ";

      my $output = { "first_fq"   => $first,
		     "first_sai"  => $first_tmp_file,
		     "second_fq"  => $second,
		     "second_sai" => $second_tmp_file};
		     
#      EASIH::JMS::submit_job($cmd, "$first_tmp_file $second_tmp_file $first $second");
      EASIH::JMS::submit_job($cmd, $output);
    }
    else {
      my $tmp_file  = EASIH::JMS::tmp_file(".sai");
      my $cmd = "$bwa aln $align_param  -f $tmp_file $reference $input ";
      my $output = { "first_fq"   => $first,
		     "first_sai"  => $tmp_file};
      EASIH::JMS::submit_job($cmd, $output);
    }
  }
  else {

    open (my $files, $input) || die "Could not open '$input': $!\n";
    while (<$files>) {
      chomp;
      my ($file1, $file2) = split("\t", $_);
      
      if ( $file1 && $file2 ) {
	my $first_tmp_file  = EASIH::JMS::tmp_file(".sai");
	my $second_tmp_file = EASIH::JMS::tmp_file(".sai");
	my $cmd = "$bwa aln $align_param  -f $first_tmp_file  $reference $file1 ;";
	$cmd   .= "$bwa aln $align_param  -f $second_tmp_file $reference $file2 ";
	my $output = { "first_fq"   => $file1,
		       "first_sai"  => $first_tmp_file,
		       "second_fq"  => $file2,
		       "second_sai" => $second_tmp_file};
		     
	EASIH::JMS::submit_job($cmd, $output);
      }
      else {
	my $tmp_file  = EASIH::JMS::tmp_file(".sai");
	my $cmd = "$bwa aln $align_param  -f $tmp_file $reference $file1 ";
	my $output = { "first_fq"   => $file1,
		       "first_sai"  => $tmp_file};
	EASIH::JMS::submit_job($cmd, $output);
      }
    }
  }
}

sub call_indels {
  my ($input) = @_;

  my @names = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$input': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }

  foreach my $name ( @names ) {

    my $tmp_file = EASIH::JMS::tmp_file(".indels");
    my $cmd = "$gatk -T IndelGenotyperV2 -R $reference -O $tmp_file -I $bam_file -L $name ";
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}		     


# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub merge_indels {
  my (@inputs) = @_;


#  print "MERGE :: @inputs \n";

  my $cmd = "cat  @inputs > $report";
  EASIH::JMS::submit_system_job("cat  @inputs > $report.indels");
  
}


# 
# Ensure that the reference and it auxiliary files are all present.
# 
# Kim Brugger (02 Aug 2010)
sub validate_input {
  
  my @errors;
  my @info;

  # Things related to the reference sequence being used.
  
  push @errors, "GATK expects references to end with 'fasta'." 
      if ( $reference !~ /fasta\z/);

  my ($dir, $basename, $postfix) = $reference =~ /^(.*)\/(.*?)\.(.*)/;
  
  push @errors, "GATK expects and references dict file (made with Picard), please see the GATK wiki\n" 
      if ( ! -e "$dir/$basename.dict");
  
  push @errors, "Platform must be either SOLEXA or SOLID not '$platform'" if ( $platform ne "SOLEXA" && $platform ne 'SOLID');

  # print the messages and die if critical ones.
  die join("\n", @errors) . "\n"   if ( @errors );
  print  join("\n", @info) . "\n"   if ( @info );
}


# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {

  $0 =~ s/.*\///;
  print "USAGE: $0 -b [am file] -R [eference genome]  -o[ut prefix] -p[latform: illumina or solid]\n";
  exit;

}
