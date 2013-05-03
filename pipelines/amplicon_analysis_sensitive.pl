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

#use lib '/home/kb468/projects/BRC_exomes/easih-toolbox/modules';
#use lib '/home/kb468/projects/BRC_exomes/easih-pipeline/modules';

use lib '/software/installed/easih-toolbox/modules';
use lib '/software/installed/easih-pipeline/modules';

use EASIH::Pipeline;
use EASIH::Pipeline::Misc;


my $opts = '1:2:b:d:De:f:hH:I:lL:mM:n:No:p:Q:Pr:R:sS:vV';
my %opts;
getopts($opts, \%opts);

usage() if ( $opts{h});

my $username = scalar getpwuid $<;

# if using standard naming, this is a lot easier.
if ( $opts{Q} ) {

  if ($opts{Q} =~ /trimmed/) {

    $opts{Q} =~ s/\.[1|2].trimmed.fq.*//;

    $opts{'1'} = join("", sort(glob("$opts{Q}*.1.trimmed.fq"), glob("$opts{Q}*.1.trimmed.fq.gz")));
    $opts{'2'} = join("", sort(glob("$opts{Q}*.2.trimmed.fq"), glob("$opts{Q}*.2.trimmed.fq.gz")));
  }
  else {
    $opts{Q} =~ s/\.[1|2].fq.*//;
    $opts{'1'} = join("", sort(glob("$opts{Q}*.1.fq"), glob("$opts{Q}*.1.fq.gz")));
    $opts{'2'} = join("", sort(glob("$opts{Q}*.2.fq"), glob("$opts{Q}*.2.fq.gz")));
  }

  $opts{'L'} = "$opts{Q}.log";
  $opts{'o'} = "$opts{Q}";
  $opts{'l'} = 1;
  $opts{'m'} = 1;

}  

#print Dumper( \%opts );

my $first          = $opts{'1'}     || usage();
#$first             = [split(",", $first)];
my $second         = $opts{'2'}     || usage();
#$second            = [split(",", $second)];
my $log            = $opts{'L'};
our $report        = $opts{'o'}     || usage();
my $platform       = 'ILLUMINA';
my $small_reference  = $opts{'r'}     || usage();
my $reference   = $opts{'R'} || '/data/refs/human_1kg/human_g1k_v37.fasta';

my $out = $opts{o} || $first;

#$out =~ s/^([A-Z]\d{6,7}).*/$1/;
$out =~ s/.1.fq//;
$out =~ s/.gz//;

my $bam_file      = "$report.bam";
my $vcf_file      = "$report.vcf";
my $host_cpus     = nr_of_cpus();

my $freeze_file = "$opts{'o'}.maris";
system "mv $freeze_file $freeze_file.backup"  if ( -e $freeze_file );

EASIH::Pipeline::freeze_file($freeze_file);



#EASIH::Pipeline::verbosity(100) if ( $opts{v});

open (*STDOUT, ">> $log") || die "Could not open '$log': $!\n" if ( $log );

# set platform specific bwa aln parameters
my $align_param .= " -q 15 ";

my $bwa             = EASIH::Pipeline::Misc::find_program('bwa_0.6.2-tpx');
my $samtools        = EASIH::Pipeline::Misc::find_program('samtools');
my $gatk            = EASIH::Pipeline::Misc::find_program('gatk_1.6.5');
my $picard          = EASIH::Pipeline::Misc::find_program('picard');


#EASIH::Pipeline::verbosity(10);
#EASIH::Pipeline::backend('Darwin');
EASIH::Pipeline::backend('Local');
EASIH::Pipeline::max_jobs( $host_cpus );
EASIH::Pipeline::max_retry(0);


EASIH::Pipeline::add_start_step('bwasw_aln');
EASIH::Pipeline::add_step('bwasw_aln', 'bam_sort');
EASIH::Pipeline::add_step('bam_sort', 'bam_index');
#EASIH::Pipeline::add_step('bam_index', 'UnifiedGenotyper');
EASIH::Pipeline::add_step('bam_index', 'amplicon_miner');
#EASIH::Pipeline::add_merge_step('UnifiedGenotyper','finished');
EASIH::Pipeline::add_merge_step('amplicon_miner','finished');

EASIH::Pipeline::print_flow();

#exit;
&EASIH::Pipeline::run();

&EASIH::Pipeline::store_state();

my $extra_report = "1 ==> $first\n";
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


EASIH::Pipeline::delete_tmp_files();

my %file2bam;



# 
# 
# 
# Kim Brugger (03 Aug 2011)
sub finished {

  system "touch $out.done";
  
}

sub bwasw_aln {
  my ($input) = @_;

  my $tmp_file = EASIH::Pipeline::tmp_file(".sam");

  my $cmd;
#  $cmd = "$bwa bwasw  -b9 -q16 -r1 -w100 $small_reference $first $second > $tmp_file 2> /dev/null";
  $cmd = "$bwa bwasw -b9 -r1 -w100 $small_reference $first $second > $tmp_file 2> /dev/null";

  print "$cmd\n";
  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_sort {
  my ($input) = @_;

  my $fixed_bam = EASIH::Pipeline::tmp_file(".bam");

  my %regions;

  open( my $sam_in,  $input ) || die "Could not open file '$input': $!\n";
  open( my $sam_out, "| $samtools view -t $reference.fai -o $fixed_bam -Sb - " ) || die "Could not open file '$fixed_bam': $!\n";
  while(<$sam_in>) {
    if (/^\@/ ) {
      print $sam_out $_ if ( !/^\@SQ/ ); 
    }
    else {
      my @F = split("\t");
      if ( $F[2] =~ /(.*?):(\d+)-(\d+)/) {
	my ($chr, $start, $end) = ($1,$2,$3);
	$regions{ $F[2] }++;
	$F[3] += $start - 1;
	$F[2] = $chr;
      }
      print $sam_out join("\t", @F);
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

  my $cmd = "$picard -T AddOrReplaceReadGroups.jar I=$fixed_bam O=$bam_file SO=coordinate CN=EASIH PL=$platform LB=$readgroup PU=1  SM=$sample VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false 2> /dev/null";

#  print "$cmd\n";

  EASIH::Pipeline::submit_job($cmd, $bam_file);
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
sub UnifiedGenotyper {
  my ($input) = @_;

  my $cmd = "$gatk -T UnifiedGenotyper -nt $host_cpus -R $reference -I $input -glm BOTH -G Standard -A AlleleBalance -stand_call_conf 30.0 -stand_emit_conf 10.0 -dcov 1000 -baq CALCULATE_AS_NECESSARY -o $vcf_file -L $out.intervals";
 

  EASIH::Pipeline::submit_job($cmd, $vcf_file);
}



# 
# 
# 
# Kim Brugger (05 Mar 2013)
sub amplicon_miner {
  my ($input) = @_;

  my $cmd = "~/projects/mosaic_miner/amplicon_miner.py $input > $out.AA.vcf";
  EASIH::Pipeline::submit_job($cmd, "$out.AA.vcf");
  
}


sub bam_rename {
  my ($input) = @_;

  EASIH::Pipeline::submit_system_job("mv $input $bam_file", $bam_file);
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
  print "USAGE: $script -1 [fastq file]  -2 [fastq file] -b[aits file (bed)] -R[eference genome] -d[bsnp rod] -o[ut prefix] -p[latform: illumina or pgm]\n";
  
  print "\nor extrapolate the standard <fq, log, out names> with the -Q flag\n";
  print "EXAMPLE: $script -Q [base name] -b[aits file (bed)] -R[eference genome] -d[bsnp rod] -p[latform: illumina or pgm]\n";

  print "\n";
  print "extra flags: -D[isable Smith-Waterman for the unmapped mate]\n";
  print "extra flags: -e[mail address, default: $username\@cam.ac.uk]\n";
  print "extra flags: -f[ilter: wgs,wgs-low,exon,exon-low. Default= exon] \n";
  print "extra flags: -I[nsert size (speeds up the run)]\n";
  print "extra flags: -L[og file, default is STDOUT]\n";  
  print "extra flags: -m[ark duplicates (always done for paired ends]\n";
  print "extra flags: -M[in depth for snps, defaults: normal=20 low=5]\n";
  print "extra flags: -N[o splitting of fastq file(s)]\n";
  print "extra flags: -n[ entries pr split-file. default: 10000000]\n";
  print "extra flags: -r[ead group]\n"; 
  print "\n";

  print "easih-pipeline: " . &EASIH::Pipeline::version() . "\n";

#  use EASIH::Toolbox;
#  print "easih-toolbox: " . &EASIH::Toolbox::version() . "\n";

  exit;

}
