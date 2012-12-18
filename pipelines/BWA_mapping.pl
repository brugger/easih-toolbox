#!/usr/bin/perl 
# 
# Map/realign/indels/snps integrated pipeline
# 
# 
# Kim Brugger (27 Jul 2010), contact: kim.brugger@easih.ac.uk

# BEGIN {
#   use vars qw/$path/; 
#   $path = $0;
#   if ($path =~ /.*\//) {
#     $path =~ s/(.*\/).*/$1/;
#   }
#   else {
#     $path = "./";
#   }
#   print "$path\n";
#   push @INC, $path;
# }

use strict;
use warnings;
use Data::Dumper;

use Getopt::Std;

use lib '/home/kb468/projects/BRC_exomes/easih-toolbox/modules';
use lib '/home/kb468/projects/BRC_exomes/easih-pipeline/modules';

use EASIH::Pipeline;
use EASIH::Pipeline::Misc;

my $opts = '1:2:d:De:f:hH:I:lL:mM:n:No:p:Q:p:Pr:R:sS:vV';
my %opts;
getopts($opts, \%opts);

usage() if ( $opts{h});
my $hard_reset    = $opts{'H'};
my $soft_reset    = $opts{'S'};

if ( $soft_reset ) {
  print "Doing a soft reset/restart\n";
  &EASIH::Pipeline::reset($soft_reset);
  getopts($opts, \%opts);
}
elsif ( $hard_reset ) {
  print "Doing a hard reset/restart\n";
  &EASIH::Pipeline::hard_reset($hard_reset);
  getopts($opts, \%opts);
}

my $username = scalar getpwuid $<;

# if using standard naming, this is a lot easier.
if ( $opts{Q} ) {

  $opts{Q} =~ s/_.*//;
  $opts{Q} =~ s/\.[1,2].fq.*//;

  $opts{'1'} = join(",", sort(glob("$opts{Q}*.1.fq"), glob("$opts{Q}*.1.fq.gz")));
  $opts{'2'} = join(",", sort(glob("$opts{Q}*.2.fq"), glob("$opts{Q}*.2.fq.gz")));

  $opts{'L'} = "$opts{Q}.log";
  $opts{'E'} = "$opts{Q}.error.log";
  $opts{'o'} = "$opts{Q}";
  $opts{'l'} = 1;
  $opts{'m'} = 1;

}  

my $first          = $opts{'1'}     || usage();
$first             = [split(",", $first)];
my $second         = $opts{'2'}     || usage();
$second            = [split(",", $second)];
my $email          = $opts{'e'}     || "$username\@cam.ac.uk";
my $loose_mapping  = $opts{'l'}     || 0;
my $log            = $opts{'L'};
my $error_log      = $opts{'E'};
my $mark_dup       = $opts{'m'};
our $report        = $opts{'o'}     || usage();
my $platform       = $opts{'p'}     || 'ILLUMINA';
my $reference      = $opts{'R'}     || usage();

my $bam_file      = "$report.bam";
my $vcf_file      = "$report.vcf";
my $host_cpus      = nr_of_cpus();

my $RUN_ON_SGE   = 1;

my $freeze_file = "$opts{'o'}.maris";
system "mv $freeze_file $freeze_file.backup"  if ( -e $freeze_file );

EASIH::Pipeline::freeze_file($freeze_file);



#EASIH::Pipeline::verbosity(100) if ( $opts{v});

open (*STDOUT, ">> $log")       || die "Could not open '$log': $!\n" if ( $log );
open (*STDERR, ">> $error_log") || die "Could not open '$error_log': $!\n" if ( $log );

# set platform specific bwa aln parameters
my $align_param .= " -q 15 ";
# and loose mapping
$align_param    .= " -e5 "     if ( $loose_mapping);

my $bwa             = EASIH::Pipeline::Misc::find_program('bwa_0.6.2-tpx');
my $samtools        = EASIH::Pipeline::Misc::find_program('samtools');
my $gatk            = EASIH::Pipeline::Misc::find_program('gatk2');
my $picard          = EASIH::Pipeline::Misc::find_program('picard');

my %run_ids         = ();

#validate_input();

#EASIH::Pipeline::verbosity(10);
#EASIH::Pipeline::backend('Darwin');
EASIH::Pipeline::backend('Local');
EASIH::Pipeline::backend('SGE') if ( $RUN_ON_SGE );
EASIH::Pipeline::max_jobs( $host_cpus ) if ( !$RUN_ON_SGE );;
EASIH::Pipeline::max_retry( 3 );

EASIH::Pipeline::add_start_step('bwa_aln');
EASIH::Pipeline::add_step('bwa_aln', 'bwa_sampe', );
EASIH::Pipeline::add_step('bwa_sampe', 'bam_sort');
EASIH::Pipeline::add_merge_step('bam_sort', 'bam_merge');
EASIH::Pipeline::add_step('bam_merge', 'bam_clean');
EASIH::Pipeline::add_step('bam_clean', 'mark_dups');
EASIH::Pipeline::add_step('mark_dups', 'bam_index');
EASIH::Pipeline::add_step('run_stats', 'finished');
EASIH::Pipeline::add_step('mark_dups', 'run_stats');
EASIH::Pipeline::add_step('bam_index', 'finished');

EASIH::Pipeline::print_flow();

EASIH::Pipeline::max_jobs( 1 ) if ( !$RUN_ON_SGE );

#exit;
&EASIH::Pipeline::run();

&EASIH::Pipeline::store_state();

my $extra_report = "1 ==> @$first\n";
$extra_report .= "bamfile ==> $bam_file\n";
$extra_report .= "snp_file ==> $report.snps\n";
$extra_report .= "indel_file ==> $report.indel\n";
$extra_report .= "easih-pipeline: " . EASIH::Pipeline::version() . "\n";

$extra_report .= "align_param ==> $align_param\n";
$extra_report .= "Binaries used..\n";
$extra_report .= "BWA: " . EASIH::Pipeline::Misc::bwa_version( $bwa ) . "\n";
$extra_report .= "Samtools: " . samtools_version() ."\n";
#$extra_report .= "GATK: " .`$gatk --version`;
$extra_report .= "Picard: " . picard_version() ."\n";
$extra_report .= "Command line: $0 ".EASIH::Pipeline::args() ."\n";


#EASIH::Pipeline::mail_report($email, $bam_file, $extra_report);


#EASIH::Pipeline::delete_tmp_files();

my %file2bam;

sub bwa_aln {
  my ($input) = @_;

  foreach my $first_file ( @$first ) {
    my $second_file = shift @$second;

    run_id( $first_file );
    
    my $tmp_sai1  = EASIH::Pipeline::tmp_file(".sai");
    my $tmp_sai2  = EASIH::Pipeline::tmp_file(".sai");
    
    my $cmd = "$bwa aln -t $host_cpus $align_param -f $tmp_sai1 $reference $first_file ; ";
    $cmd   .= "$bwa aln -t $host_cpus $align_param -f $tmp_sai2 $reference $second_file";
    
    my $output = { "first_fq"   => $first_file,
		   "first_sai"  => $tmp_sai1,
		   "second_fq"  => $second_file,
		   "second_sai" => $tmp_sai2};
    
    EASIH::Pipeline::submit_job($cmd, $output, "t=$host_cpus");
  }

}

# 
# 
# 
# Kim Brugger (10 Dec 2012)
sub run_id {
  my ( $fq_file ) = @_;

  my $header_line = `zcat -f $fq_file | head -n 1`;
  chomp ( $header_line );
  $header_line =~ s/^\@//;
  my ($run_id, $lane,  @rest) = split(/:/, $header_line);

  $run_ids{ $fq_file } = "$run_id.$lane";
}


sub bwa_sampe {
  my ($input) = @_;

  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");

  my $cmd;
  $cmd = "$bwa sampe -t $host_cpus -P $reference $$input{first_sai} $$input{second_sai} $$input{first_fq} $$input{second_fq} -f $tmp_file";

  $file2bam{ $tmp_file} = $$input{ 'first_fq' };

  EASIH::Pipeline::submit_job($cmd, $tmp_file, "t=$host_cpus");
}



# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_sort {
  my ($input) = @_;

  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");

  my $readgroup =   $file2bam{ $input };
  $readgroup =~ s/\.gz//;
  $readgroup =~ s/\.[1|2].fq//;
  $readgroup =~ s/\_\d+//;
  
  my $sample = $readgroup;
  $sample =~ s/\..*//;
  $sample =~ s/_\d*//;

  my $run_id    = $run_ids{ $file2bam{ $input } };


  my $cmd = "$picard -T AddOrReplaceReadGroups.jar I=$input O=$tmp_file SORT_ORDER=coordinate CN=EASIH PL=$platform LB=$readgroup PU=$run_id  SM=$sample VALIDATION_STRINGENCY=SILENT ";

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}



# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_merge {
  my (@inputs) = @_;

  EASIH::Pipeline::max_jobs( $host_cpus ) if ( ! $RUN_ON_SGE );


  @inputs = @{$inputs[0]} if ( @inputs == 1 && ref($inputs[0]) eq "ARRAY" );


  my $tmp_file = EASIH::Pipeline::tmp_file(".merged.bam");

  if (@inputs == 1 ) {
    EASIH::Pipeline::submit_system_job("mv @inputs $tmp_file", $tmp_file);
  }
  else {

    # remove empty files as they crash the merging step.
    my @non_empty_files;
    foreach my $input ( @inputs ) {
      push @non_empty_files, $input if ( ! -z $input );
    }

    my $username = scalar getpwuid $<;

    my $cmd = "$picard -T MergeSamFiles USE_THREADING=true O=$tmp_file  I= " . join(" I= ", @non_empty_files) . " VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/ CREATE_INDEX=false";
    EASIH::Pipeline::submit_job($cmd, $tmp_file);
  }
  
}


# 
# 
# 
# Kim Brugger (26 Jun 2012)
sub bam_index {
  my ($input) = @_;

  my $cmd = "$samtools index $input";
  EASIH::Pipeline::submit_job($cmd, $input);
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub mark_dups {
  my ($input) = @_;

  my $username = scalar getpwuid $<;
  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");
  my $metrix_file = EASIH::Pipeline::tmp_file(".mtx");
  my $cmd = "$picard -T MarkDuplicates  I= $input O= $bam_file  M= $metrix_file VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/ MAX_RECORDS_IN_RAM=500000";
  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}



# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_clean {
  my ($input) = @_;

  my $username = scalar getpwuid $<;
  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");
  my $cmd = "$picard -T CleanSam I= $input O= $tmp_file VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/ MAX_RECORDS_IN_RAM=500000";
  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}


sub run_stats {
  my ($input) = @_;
  
  EASIH::Pipeline::submit_job("$samtools flagstat $report.bam > $report.flagstat");
  EASIH::Pipeline::submit_job("md5sum $bam_file > $bam_file.md5 ");
  EASIH::Pipeline::submit_job("md5sum $report.bam.bai > $report.bam.bai.md5 ");
  EASIH::Pipeline::submit_job("md5sum $report.maris > $report.maris.md5");
}


sub rename {
  my ($input) = @_;

  EASIH::Pipeline::submit_system_job("mv $input $bam_file", $bam_file);
}


# 
# 
# 
# Kim Brugger (03 Aug 2011)
sub finished {

  open(my $out, "> $report.done") || die "Could not write to '$report.done': $!\n";
  print $out "$bam_file\n";
  print $out "$bam_file.md5\n";
  print $out "$bam_file.bai\n";
  print $out "$bam_file.bai.md5\n";
  print $out "$report.maris\n";
  print $out "$report.maris.md5\n";
  print $out "$report.flagstat\n";
  close( $out);
  
}




# 
# Ensure that the reference and it auxiliary files are all present.
# 
# Kim Brugger (02 Aug 2010)
sub validate_input {
  
  my @errors;
  my @info;

  
  foreach my $first_file ( @$first ) {
    push @errors, "$first_file does not exist"  if ( ! -e $first_file );
  }

  foreach my $second_file ( @$second ) {
    push @errors, "$second_file does not exist"  if ( ! $second_file );
  }

  # Things related to the reference sequence being used.
  
  push @errors, "GATK expects references to end with 'fasta'." 
      if ( $reference !~ /fasta\z/);

  my ($dir, $basename, $postfix) = (".","","");
  if ( $reference =~ /\//) {
    ($dir, $basename, $postfix) = $reference =~ /^(.*)\/(.*?)\.(.*)/;
  }
  else {
    ($basename, $postfix) = $reference =~ /^(.*?)\.(.*)/;
  }
  
  push @errors, "GATK expects and references dict file (made with Picard), please see the GATK wiki $dir/$basename.dict\n" 
      if ( ! -e "$dir/$basename.dict");
  
  my @bwa_postfixes = ('amb', 'ann', 'bwt', 'fai','pac', 'rbwt', 'rpac', 'rsa', 'sa');

 
  foreach my $bwa_postfix ( @bwa_postfixes ) {
    push @errors, "$reference.$bwa_postfix does not exists. Did you run bwa index on $reference?"
	if ( ! -e "$reference.$bwa_postfix");
  }


  # Check that the bam_file ends with bam, or add it
  if ( $bam_file !~ /bam\z/) {
    push @info, "Added bam postfix so '$bam_file' becomes '$bam_file.bam'";
    $bam_file .= ".bam";
  }

  push @errors, "Platform must be either ILLUMINA or PGM not '$platform'" if ( $platform ne "ILLUMINA" && $platform ne 'PGM');


  # print the messages and die if critical ones.
  die join("\n", @errors) . "\n"   if ( @errors );
  print  join("\n", @info) . "\n"   if ( @info );
}



# 
# 
# 
# Kim Brugger (13 Jan 2011)
sub nr_of_cpus {

  my $cpus = `cat /proc/cpuinfo | egrep ^proc | wc -l`;
  chomp $cpus;
  return $cpus;
}



sub samtools_version {
  open( my $pipe, "$samtools 2>&1 | head -n 3 |") || die "Could not open samtools pipe: $!\n";
  <$pipe>;
  <$pipe>;
  my $version = <$pipe>;
  chomp( $version);

  return( $version);
}



# 
# 
# 
# Kim Brugger (08 Nov 2010)
sub picard_version {

  my $version = `$picard -v`;
  chomp($version);
  return $version;
}


# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {

  my $script = $0;
  $script =~ s/.*\///;
  print "USAGE: $script -1 [fastq file]  -2 [fastq file] -l[oose mapping] -r[eference directory] -o[ut prefix]\n";
  
  print "\nor extrapolate the standard <fq, log, out names> with the -Q flag\n";
  print "EXAMPLE: $script -Q [base name]  -r[eference directory] \n";

  print "\n";
  print "extra flags: -D[isable Smith-Waterman for the unmapped mate]\n";
  print "extra flags: -e[mail address, default: $username\@cam.ac.uk]\n";
  print "extra flags: -H[ard reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "extra flags: -L[og file, default is STDOUT]\n";  
  print "extra flags: -S[oft reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "\n";

  print "easih-pipeline: " . &EASIH::Pipeline::version() . "\n";

  use EASIH::Toolbox;
  print "easih-toolbox: " . &EASIH::Toolbox::version() . "\n";

  exit;

}
