#!/usr/bin/perl 
# 
# pipeline for few genes analysis, mainly for clinical use.
# 
# 
# Kim Brugger (17 Jan 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

BEGIN {
  use vars qw/$path/; 
  $path = $0;
  if ($path =~ /.*\//) {
    $path =~ s/(.*\/).*/$1/;
  }
  else {
    $path = "./";
  }
  push @INC, $path;
}

use Getopt::Std;

#use lib '/home/cjp64/git/easih-pipeline/modules';
use lib '/home/kb468/projects/BRC_exomes/easih-pipeline/modules';
use lib '/home/kb468/projects/BRC_exomes/easih-toolbox/modules';

use EASIH::Misc;
use EASIH::Pipeline;

my $opts = '1:2:R:r:p:B:o:d:';
my %opts;
getopts($opts, \%opts);

#usage() if ( $opts{h});


my $first         = $opts{'1'}     || usage();
my $second         = $opts{'2'}     || usage();
my $platform = 'ILLUMINA';
my $reference        = $opts{'R'}     || usage();
my $small_reference  = $opts{'r'}     || usage();
my $leeway           = $opts{'l'}     || 50;
my $dbsnp            = $opts{d} || usage();
my $baits            = $opts{B}  || usage();


my $samtools    = EASIH::Misc::find_program('samtools');
my $smalt       = EASIH::Misc::find_program('smalt');
my $picard      = EASIH::Misc::find_program('picard');
my $gatk        = EASIH::Misc::find_program('gatk_1.3-14');

# set platform specific bwa aln parameters
my $align_param .= " -q 15 -e 5";

my $bwa             = EASIH::Misc::find_program('bwa_0.6.1-tpx');


my $smalt_ref = $reference;
$smalt_ref =~ s/^(.*)\..*/$1/;

my $out = $opts{o} || $first;

#$out =~ s/^([A-Z]\d{6,7}).*/$1/;
$out =~ s/.1.fq//;
$out =~ s/.gz//;

my $bam_file = "$out.bam";
my $vcf_file = "$out.vcf";

my %regions;

my $username = scalar getpwuid $<;
my $host_cpus      = nr_of_cpus();
#$host_cpus = 4;

EASIH::Pipeline::add_start_step('bwa_aln');
EASIH::Pipeline::add_step('bwa_aln', 'bwa_sampe', );
EASIH::Pipeline::add_step('bwa_sampe', 'bam_sort');
EASIH::Pipeline::add_step('bam_sort', 'mark_dups');
EASIH::Pipeline::add_step('mark_dups', 'bam_index');
EASIH::Pipeline::add_step('bam_index', 'realign_targets');
EASIH::Pipeline::add_step('realign_targets', 'bam_realign');
EASIH::Pipeline::add_merge_step('bam_realign', 'bam_merge2' );
EASIH::Pipeline::add_step('bam_merge2', 'bam_index2');
EASIH::Pipeline::add_step('bam_index2', 'count_covariates');
EASIH::Pipeline::add_step('count_covariates', 'table_recalibration');
EASIH::Pipeline::add_step('table_recalibration', 'bam_index3');

EASIH::Pipeline::add_step('bam_index3', 'UnifiedGenotyper');

EASIH::Pipeline::add_step('UnifiedGenotyper', 'annotate_VCF');
EASIH::Pipeline::add_step('annotate_VCF', 'vcf2xls');

EASIH::Pipeline::backend('Local');
EASIH::Pipeline::max_jobs( $host_cpus );
EASIH::Pipeline::max_retry(0);

EASIH::Pipeline::print_flow();
EASIH::Pipeline::run();

my %file2bam;


sub bwa_aln {
  my ($input) = @_;

  my $tmp_sai1  = EASIH::Pipeline::tmp_file(".sai");
  my $tmp_sai2  = EASIH::Pipeline::tmp_file(".sai");
    
  my $cmd = "$bwa aln -t $host_cpus $align_param -f $tmp_sai1 $small_reference $first ; ";
  $cmd   .= "$bwa aln -t $host_cpus $align_param -f $tmp_sai2 $small_reference $second";
  
  my $output = { "first_fq"   => $first,
		 "first_sai"  => $tmp_sai1,
		 "second_fq"  => $second,
		 "second_sai" => $tmp_sai2};
  
  EASIH::Pipeline::submit_job($cmd, $output);

}


sub bwa_sampe {
  my ($input) = @_;

  my $tmp_file = EASIH::Pipeline::tmp_file(".sam");

  my $cmd;
  $cmd = "$bwa sampe -t $host_cpus -P $small_reference $$input{first_sai} $$input{second_sai} $$input{first_fq} $$input{second_fq} > $tmp_file";

  $file2bam{ $tmp_file} = $$input{ 'first_fq' };

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}

# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_sort {
  my ($input) = @_;

  my $fixed_bam = EASIH::Pipeline::tmp_file(".bam");
  open( my $sam_in,  $input ) || die "Could not open file '$input': $!\n";
  open( my $sam_out, "| $samtools view -t $reference.fai -Sb - -o $fixed_bam" ) || die "Could not open file '$fixed_bam': $!\n";
  while(<$sam_in>) {
    if (/^\@/ && !/^\@SQ/) {
      print $sam_out $_;
    }
    else {
      my @F = split("\t");
      if ( $F[2] =~ /(.*?):(\d+)-(\d+)/) {
	my ($chr, $start, $end) = ($1,$2,$3);
	$regions{ $F[2] }++;
	$F[3] += $start - 1;
	$F[2] = $chr;
	print $sam_out join("\t", @F);
      }
    }
  }
  close($sam_out);
  close($sam_in);

  open (my $var_out, " > $out.intervals") || die "Could not open file '$out.intervals': $!\n";
  foreach my $region (keys %regions) {
    print $var_out "$region\n";
  }
  close( $var_out );


  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");

  my $readgroup =   $first;
  $readgroup =~ s/\.gz//;
  $readgroup =~ s/\.[1|2].fq//;
  
  my $sample = $readgroup;
  $sample =~ s/\..*//;
  $sample =~ s/_\d*//;

  my $cmd = "$picard -T AddOrReplaceReadGroups.jar I=$fixed_bam O=$tmp_file SO=coordinate CN=EASIH PL=$platform LB=$readgroup PU=1  SM=$sample VALIDATION_STRINGENCY=SILENT ";

  print "$cmd\n";

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
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
  my $cmd = "$picard -T MarkDuplicates  I= $input O= $tmp_file  M= $metrix_file VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/ MAX_RECORDS_IN_RAM=500000";
  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub count_covariates {
  my ($input) = @_;
  
  my $tmp_file = EASIH::Pipeline::tmp_file("_recal.csv.");

  my $cmd = "$gatk -T CountCovariates -R $reference -I $input -cov ReadGroupCovariate  -cov QualityScoreCovariate  -cov CycleCovariate  -cov DinucCovariate  -recalFile $tmp_file -L $out.intervals";
  $cmd .= " -knownSites $dbsnp ";

  EASIH::Pipeline::submit_job($cmd, "$input $tmp_file");
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub table_recalibration {
  my ($input) = @_;
  my ($tmp_bam, $recal) = split(" ", $input);

  my $cmd = "$gatk -T TableRecalibration  -R $reference -I $tmp_bam -recalFile $recal -baq RECALCULATE -o $bam_file -L $out.intervals"; 

  EASIH::Pipeline::submit_job($cmd, $bam_file);
}

# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub UnifiedGenotyper {
  my ($input) = @_;

  my $tmp_file = EASIH::Pipeline::tmp_file(".vcf");
  my $cmd = "$gatk -T UnifiedGenotyper -nt $host_cpus -R $reference -I $input -glm BOTH -G Standard -A AlleleBalance -stand_call_conf 30.0 -stand_emit_conf 10.0 -dcov 1000 -baq CALCULATE_AS_NECESSARY -o $vcf_file";
  $cmd .= " -L $out.intervals";
  
  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}



# 
# 
# 
# Kim Brugger (05 Jul 2012)
sub annotate_VCF {


  my $cmd = "~/projects/perl_xls/clinical_report.pl -v $vcf_file -O $out.var.csv -o $out.var_full.csv -b $bam_file -B $baits\n";
  EASIH::Pipeline::submit_job($cmd);

}



# 
# 
# 
# Kim Brugger (05 Jul 2012)
sub vcf2xls {
  my $cmd = "~/projects/perl_xls/v2.2/clinical_report_builder_2.21.pl $out.var_full.csv \n";
  EASIH::Pipeline::submit_job($cmd);
  
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
# Kim Brugger (26 Jun 2012)
sub bam_index2 {
  my ($input) = @_;

  my $cmd = "$samtools index $input";
  EASIH::Pipeline::submit_job($cmd, $input);
}



# 
# 
# 
# Kim Brugger (26 Jun 2012)
sub bam_index3 {
  my ($input) = @_;

  my $cmd = "$samtools index $input";
  EASIH::Pipeline::submit_job($cmd, $input);
}


# 
# 
# 
# Kim Brugger (26 Jun 2012)
sub bam_index4 {
  my ($input) = @_;

  my $cmd = "$samtools index $input";
  EASIH::Pipeline::submit_job($cmd, $input);
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub realign_targets {
  my ($input) = @_;

  my $tmp_file = EASIH::Pipeline::tmp_file(".intervals");
  my $cmd = "$gatk -T RealignerTargetCreator -R $reference -o $tmp_file -L $out.intervals -I $input";
  EASIH::Pipeline::submit_job($cmd, "$tmp_file $input");
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_realign {
  my ($input) = @_;

  my ($interval_file, $tmp_bam_file) = split(" ", $input);

  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");
  my $cmd;
  $cmd = "$gatk -T IndelRealigner -targetIntervals $interval_file -L $out.intervals -o $tmp_file -R $reference -I $tmp_bam_file ";

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_merge2 {
  my (@inputs) = @_;


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

    my $cmd = "$picard -T MergeSamFiles USE_THREADING=true O=$tmp_file  I= " . join(" I= ", @non_empty_files) . " VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/";
    EASIH::Pipeline::submit_job($cmd, $tmp_file);
  }
  
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

# 
# 
# 
# Kim Brugger (17 Jan 2012)
sub usage {


  my $script = $0;
  $script =~ s/.*\///;
  print "USAGE: $script -1 [fastq file]  -2 [fastq file]  -r[eference directory] -o[ut prefix]\n";
  
  print "\nor extrapolate the standard <fq, log, out names> with the -Q flag\n";
  print "EXAMPLE: $script -Q [base name] -l[oose mapping] -R[eference genome] -d[bsnp rod] -p[latform: illumina or pgm]\n";

  print "\n";
  print "extra flags: -D[isable Smith-Waterman for the unmapped mate]\n";
  print "extra flags: -e[mail address, default: $username\@cam.ac.uk]\n";
  print "extra flags: -f[ilter: wgs,wgs-low,exon,exon-low. Default= exon] \n";
  print "extra flags: -H[ard reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "extra flags: -I[nsert size (speeds up the run)]\n";
  print "extra flags: -L[og file, default is STDOUT]\n";  
  print "extra flags: -m[ark duplicates (always done for paired ends]\n";
  print "extra flags: -M[in depth for snps, defaults: normal=20 low=5]\n";
  print "extra flags: -N[o splitting of fastq file(s)]\n";
  print "extra flags: -n[ entries pr split-file. default: 10000000]\n";
  print "extra flags: -r[ead group]\n"; 
  print "extra flags: -S[oft reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "\n";

  print "easih-pipeline: " . &EASIH::Pipeline::version() . "\n";

  use EASIH::Toolbox;
  print "easih-toolbox: " . &EASIH::Toolbox::version() . "\n";

  exit;

}




