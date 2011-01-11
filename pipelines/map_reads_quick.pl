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



our %analysis = ('fastq-split'      => { function   => 'fastq_split',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=02:00:00"},
		 
		 'std-aln'          => { function   => 'bwa_aln',
					 func_param => ' ',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=12:00:00"},
		 
		 'std-generate'      => { function   => 'bwa_generate',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},

		 'std-tag_sam'       => { function   => 'sam_add_tags',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},
		 
		 'std-sam2bam'       => { function   => 'EASIH::JMS::Samtools::sam2bam',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=1000mb,walltime=10:00:00"},
		 
		 'std-merge'         => { function   => 'EASIH::JMS::Picard::merge',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",
					 sync       => 1},

		 'std-sort'    => { function   => 'EASIH::JMS::Picard::sort',
				    hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=20000mb,walltime=08:00:00"},
		 
		 
		 'std-rename'     => { function   => 'rename'},

		 
		 'std-index'   =>{ function   => 'EASIH::JMS::Samtools::index',
				   hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
		 

		 );
		     


our %flow = ( 'csfasta2fastq'    => 'std-aln',
	      'fastq-split'      => 'std-aln',

	      'std-aln'          => 'std-generate',
	      'std-generate'     => 'std-tag_sam',
	      'std-tag_sam'      => 'std-sam2bam',
	      'std-sam2bam'      => 'std-merge',
	      'std-merge'        => 'std-sort',
	      'std-sort'         => 'std-rename',
	      'std-rename'       => 'std-index'
	      );

#EASIH::JMS::no_store();
#EASIH::JMS::print_flow('fastq-split');

my %opts;
getopts('1:2:nm:R:o:r:p:hH:lL:S:', \%opts);


usage() if ( $opts{h});
my $hard_reset    = $opts{'H'};
my $soft_reset    = $opts{'S'};

if ( $soft_reset ) {
  &EASIH::JMS::reset($soft_reset);
  getopts('1:2:nm:R:o:r:p:hH:lL:S:', \%opts);
}
elsif ( $hard_reset ) {
  &EASIH::JMS::hard_reset($hard_reset);
  getopts('1:2:nm:R:o:r:p:hH:lL:S:', \%opts);
}


my $first       = $opts{'1'} || usage();
my $second      = $opts{'2'};
my $no_split    = $opts{'n'} || 0;
my $split       = $opts{'m'} || 5000000;

my $log         = $opts{'L'};
open (*STDOUT, ">> $log") || die "Could not open '$log': $!\n" if ( $log );

my $reference   = $opts{'R'} || usage();
my $align_param = ' ';
my $report      = $opts{'o'} || usage();
my $bam_file    = "$report.bam" || usage();
my $loose_mapping = $opts{'l'}    || 0;

my $readgroup   = $opts{'r'} || $report;
my $platform    = uc($opts{'p'}) || usage();
$platform = 'SOLEXA'      if ( $platform eq 'ILLUMINA');

# set platform specific bwa aln parameters
$align_param .= " -c "      if ( $platform eq "SOLID");
$align_param .= " -q 15 "   if ( $platform eq "SOLEXA");
$align_param .= " -e5 -t5 " if ( $loose_mapping);

my $bwa          = EASIH::JMS::Misc::find_program('bwa');
my $fq_split     = EASIH::JMS::Misc::find_program('fastq_split.pl');
my $samtools     = EASIH::JMS::Misc::find_program('samtools');
my $tag_sam      = EASIH::JMS::Misc::find_program('tag_sam.pl');
my $sam2fq       = EASIH::JMS::Misc::find_program('sam2fastq.pl');

validate_input();

#EASIH::JMS::verbosity(10);
EASIH::JMS::backend('Darwin');
#EASIH::JMS::backend('Kluster');
EASIH::JMS::max_retry(0);


if ( $no_split ) {
  &EASIH::JMS::run('std-aln');
}
else {
  &EASIH::JMS::run('fastq-split');
}

&EASIH::JMS::store_state();

my $extra_report = "1 ==> $first\n";
$extra_report .= "2 ==> $second\n" if ( $second );
$extra_report .= "bamfile ==> $bam_file\n";
$extra_report .= "snp_file ==> $report.snps\n";
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
      my $cmd = "$bwa aln $align_param  -f $tmp_file $reference $first ";
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


sub bwa_generate {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".sam");

  my $cmd;
  if (defined($$input{'second_sai'}) ) {
    $cmd = "$bwa sampe -f $tmp_file $reference $$input{first_sai} $$input{second_sai} $$input{first_fq} $$input{second_fq}";
  
  }
  else {
    $cmd = "$bwa samse -f $tmp_file $reference $$input{first_sai} $$input{first_fq}";
  }

  EASIH::JMS::submit_job($cmd, $tmp_file);
}

# 
# 
# 
# Kim Brugger (04 Aug 2010)
sub bam2fq {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".fq");
  my $list_tmp_file = EASIH::JMS::tmp_file(".list");
  my $cmd = "$samtools view $input | $sam2fq -2 $tmp_file > $list_tmp_file";  
  EASIH::JMS::submit_job($cmd,  $list_tmp_file);
}



sub bwa_aln_loose {
  my ($input) = @_;

  if ( $no_split ) {
    
    if ( $first && $second ) {
      my $first_tmp_file  = EASIH::JMS::tmp_file(".sai");
      my $second_tmp_file = EASIH::JMS::tmp_file(".sai");
      my $cmd = "$bwa aln $align_param -e5 -t5 -f $first_tmp_file  $reference $first ;";
      $cmd   .= "$bwa aln $align_param -e5 -t5 -f $second_tmp_file $reference $second ";

      my $output = { "first_fq"   => $first,
		     "first_sai"  => $first_tmp_file,
		     "second_fq"  => $second,
		     "second_sai" => $second_tmp_file};
		     
      EASIH::JMS::submit_job($cmd, $output);
    }
    else {
      my $tmp_file  = EASIH::JMS::tmp_file(".sai");
      my $cmd = "$bwa aln $align_param -e5 -t5 -f $tmp_file $reference $input ";
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
	my $cmd = "$bwa aln $align_param -e5 -t5 -f $first_tmp_file  $reference $file1 ;";
	$cmd   .= "$bwa aln $align_param -e5 -t5 -f $second_tmp_file $reference $file2 ";
	my $output = { "first_fq"   => $file1,
		       "first_sai"  => $first_tmp_file,
		       "second_fq"  => $file2,
		       "second_sai" => $second_tmp_file};
		     
	EASIH::JMS::submit_job($cmd, $output);
      }
      else {
	my $tmp_file  = EASIH::JMS::tmp_file(".sai");
	my $cmd = "$bwa aln $align_param -e5 -t5 -f $tmp_file $reference $file1 ";
	my $output = { "first_fq"   => $file1,
		       "first_sai"  => $tmp_file};
	EASIH::JMS::submit_job($cmd, $output);
      }
    }
  }

}



# 
# Adds readgroup & aligner tags to the sam-file.
# 
# Kim Brugger (22 Jul 2010)
sub sam_add_tags {
  my ($input) = @_;

  if ( ! $readgroup ) {
    $readgroup = $first;
    $readgroup =~ s/.fastq//;
    $readgroup =~ s/.fq//;
    $readgroup =~ s/.gz//;
  }

  my $cmd = "$tag_sam -R $input -p $platform -r $readgroup -a bwa -A '$align_param' ";
  EASIH::JMS::submit_job($cmd, $input);
}




sub rename {
  my ($input) = @_;

  EASIH::JMS::submit_system_job("mv $input $bam_file", $bam_file);
}





# 
# Ensure that the reference and it auxiliary files are all present.
# 
# Kim Brugger (02 Aug 2010)
sub validate_input {
  
  my @errors;
  my @info;

  # Things related to the reference sequence being used.
  

  my ($dir, $basename, $postfix) = (".","","");
  if ( $reference =~ /\//) {
    ($dir, $basename, $postfix) = $reference =~ /^(.*)\/(.*?)\.(.*)/;
  }
  else {
    ($basename, $postfix) = $reference =~ /^(.*?)\.(.*)/;
  }

  
  
  my @bwa_postfixes = ('amb', 'ann', 'bwt', 'fai','pac', 'rbwt', 'rpac', 'rsa', 'sa');

  push @bwa_postfixes, ( 'nt.amb', 'nt.ann', 'nt.pac')  if ( $platform eq "SOLID");
  
  foreach my $bwa_postfix ( @bwa_postfixes ) {
    push @errors, "$reference.$bwa_postfix does not exists. Did you run bwa index on $reference?"
	if ( ! -e "$reference.$bwa_postfix");
  }


  # Check that the bam_file ends with bam, or add it
  if ( $bam_file !~ /bam\z/) {
    push @info, "Added bam postfix so '$bam_file' becomes '$bam_file.bam'";
    $bam_file .= ".bam";
  }


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
  print "USAGE: $0 -1 [fastq file]  -2 [fastq file] -l[oose mapping] -n[o splitting of fastq file(s)] -R [eference genome]  -o[ut prefix] -p[latform: illumina or solid]\n";

  print "\n";
  print "extra flags: -H[ard reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "extra flags: -L[og file, default is STDOUT]\n";  
  print "extra flags: -S[oft reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "\n";
  print "easih-pipeline: " . &EASIH::JMS::version() . "\n";
  exit;

}
