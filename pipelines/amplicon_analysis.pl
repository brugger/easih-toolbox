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


my $opts = '1:2:b:d:De:f:hH:I:lL:mM:n:No:p:Q:Pr:R:sS:vV';
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

  $opts{'1'} = join(",", sort(glob("$opts{Q}*.1.fq"), glob("$opts{Q}*.1.fq.gz")));
  $opts{'2'} = join(",", sort(glob("$opts{Q}*.2.fq"), glob("$opts{Q}*.2.fq.gz")));

  $opts{'L'} = "$opts{Q}.log";
  $opts{'o'} = "$opts{Q}";
  $opts{'l'} = 1;
  $opts{'m'} = 1;

}  

my $baits          = $opts{'b'}     || usage();
my $first          = $opts{'1'}     || usage();
$first             = [split(",", $first)];
my $second         = $opts{'2'}     || usage();
$second            = [split(",", $second)];
my $email          = $opts{'e'}     || "$username\@cam.ac.uk";
my $loose_mapping  = $opts{'l'}     || 0;
my $log            = $opts{'L'};
my $mark_dup       = $opts{'m'};
our $report        = $opts{'o'}     || usage();
my $platform       = 'ILLUMINA';
my $reference_dir = $opts{'r'}     || usage();

my $reference   = glob("$reference_dir/*.fasta");
my ($dbsnp)     = grep(!/excluding/, glob("$reference_dir/dbsnp_*.vcf"));
my ($hapmap)    = glob("$reference_dir/hapmap_*.sites.vcf");
my ($omni)      = glob("$reference_dir/*omni*.sites.vcf");

my $bam_file      = "$report.bam";
my $vcf_file      = "$report.vcf";
my $host_cpus     = 2;#nr_of_cpus();

my $freeze_file = "$opts{'o'}.maris";
system "mv $freeze_file $freeze_file.backup"  if ( -e $freeze_file );

EASIH::Pipeline::freeze_file($freeze_file);


my $run_id = "MGILLUMINA4_75";

#EASIH::Pipeline::verbosity(100) if ( $opts{v});

open (*STDOUT, ">> $log") || die "Could not open '$log': $!\n" if ( $log );

# set platform specific bwa aln parameters
my $align_param .= " -q 15 ";
# and loose mapping
$align_param    .= " -e5 "     if ( $loose_mapping);

my $bwa             = EASIH::Pipeline::Misc::find_program('bwa_0.6.1-tpx');
my $samtools        = EASIH::Pipeline::Misc::find_program('samtools');
my $gatk            = EASIH::Pipeline::Misc::find_program('gatk_1.6.5');
my $picard          = EASIH::Pipeline::Misc::find_program('picard');

my $mpi_q = "/home/kb468/easih-pipeline/tools/mpiexec_queue.pl";

#validate_input();


#EASIH::Pipeline::verbosity(10);
#EASIH::Pipeline::backend('Darwin');
EASIH::Pipeline::backend('Local');
EASIH::Pipeline::max_jobs( $host_cpus );
EASIH::Pipeline::max_retry(0);

EASIH::Pipeline::add_start_step('bwa_aln');
EASIH::Pipeline::add_step('bwa_aln', 'bwa_sampe', );
EASIH::Pipeline::add_step('bwa_sampe', 'bam_sort');
EASIH::Pipeline::add_step('bam_sort', 'bam_rename');
EASIH::Pipeline::add_step('bam_rename', 'bam_index');
EASIH::Pipeline::add_step('bam_index', 'UnifiedGenotyper');

EASIH::Pipeline::print_flow();

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


EASIH::Pipeline::delete_tmp_files();

my %file2bam;

sub bwa_aln {
  my ($input) = @_;

  foreach my $first_file ( @$first ) {
    my $second_file = shift @$second;
    
    my $tmp_sai1  = EASIH::Pipeline::tmp_file(".sai");
    my $tmp_sai2  = EASIH::Pipeline::tmp_file(".sai");
    
    my $cmd = "$bwa aln -t $host_cpus $align_param -f $tmp_sai1 $reference $first_file ; ";
    $cmd   .= "$bwa aln -t $host_cpus $align_param -f $tmp_sai2 $reference $second_file";
    
    my $output = { "first_fq"   => $first_file,
		   "first_sai"  => $tmp_sai1,
		   "second_fq"  => $second_file,
		   "second_sai" => $tmp_sai2};
    
    EASIH::Pipeline::submit_job($cmd, $output);
  }

}


sub bwa_sampe {
  my ($input) = @_;

  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");

  my $cmd;
  $cmd = "$bwa sampe -t $host_cpus -P $reference $$input{first_sai} $$input{second_sai} $$input{first_fq} $$input{second_fq} | samtools view  -t $reference -Sb - > $tmp_file";

  $file2bam{ $tmp_file} = $$input{ 'first_fq' };

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
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
  
  my $sample = $readgroup;
  $sample =~ s/\..*//;
  $sample =~ s/_\d*//;

  my $cmd = "$picard -T AddOrReplaceReadGroups.jar I=$input O=$tmp_file SORT_ORDER=coordinate CN=EASIH PL=$platform LB=$readgroup PU=$run_id  SM=$sample VALIDATION_STRINGENCY=SILENT ";

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
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

  my $tmp_file = EASIH::Pipeline::tmp_file(".vcf");
  my $cmd = "$gatk -T UnifiedGenotyper -nt $host_cpus -R $reference -I $input -glm BOTH -G Standard -A AlleleBalance -stand_call_conf 30.0 -stand_emit_conf 10.0 -dcov 1000 -baq CALCULATE_AS_NECESSARY -o $vcf_file -L $baits";
  $cmd .= " --dbsnp $dbsnp";
  

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
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

  push @errors, "'$dbsnp' does not exists\n" if (! -e $dbsnp);
  push @errors, "'$dbsnp' don't end with .vcf as expected\n" if ($dbsnp !~ /.vcf\z/);

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
