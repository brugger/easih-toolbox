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


my $opts = '1:2:d:De:f:hH:I:lL:mM:n:No:p:Q:Pr:R:sS:vV';
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

  $opts{Q} =~ s/(.*?)\..*/$1/;
  

  $opts{'1'} = join(",", sort(glob("$opts{Q}*.1.fq"), glob("$opts{Q}*.1.fq.gz")));
  $opts{'2'} = join(",", sort(glob("$opts{Q}*.2.fq"), glob("$opts{Q}*.2.fq.gz")));

  $opts{'L'} = "$opts{Q}.log";
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
my $mark_dup       = $opts{'m'};
our $report        = $opts{'o'}     || usage();
my $platform       = 'ILLUMINA';

my $reference   = "/data/refs/HIV/K03455.fasta";
my $bam_file      = "$report.bam";
my $host_cpus      = nr_of_cpus();

my $freeze_file = "$opts{'o'}.maris";
system "mv $freeze_file $freeze_file.backup"  if ( -e $freeze_file );

EASIH::Pipeline::freeze_file($freeze_file);

my $run_id = "---";


open (*STDOUT, ">> $log") || die "Could not open '$log': $!\n" if ( $log );


my $samtools        = EASIH::Pipeline::Misc::find_program('samtools');
my $picard          = EASIH::Pipeline::Misc::find_program('picard');

#validate_input();


#EASIH::Pipeline::verbosity(10);
#EASIH::Pipeline::backend('Darwin');
EASIH::Pipeline::backend('Local');
EASIH::Pipeline::max_jobs( $host_cpus );
EASIH::Pipeline::max_retry(0);

EASIH::Pipeline::add_start_step('smalt');
EASIH::Pipeline::add_step('smalt', 'bam_sort');
EASIH::Pipeline::add_step('bam_sort', 'mark_dups');


EASIH::Pipeline::print_flow();

#exit;
&EASIH::Pipeline::run();

&EASIH::Pipeline::store_state();

my $extra_report = "1 ==> @$first\n";
$extra_report .= "bamfile ==> $bam_file\n";
$extra_report .= "snp_file ==> $report.snps\n";
$extra_report .= "indel_file ==> $report.indel\n";
$extra_report .= "easih-pipeline: " . EASIH::Pipeline::version() . "\n";

$extra_report .= "Binaries used..\n";
$extra_report .= "Samtools: " . samtools_version() ."\n";
#$extra_report .= "GATK: " .`$gatk --version`;
$extra_report .= "Picard: " . picard_version() ."\n";
$extra_report .= "Command line: $0 ".EASIH::Pipeline::args() ."\n";


#EASIH::Pipeline::mail_report($email, $bam_file, $extra_report);


#EASIH::Pipeline::delete_tmp_files();

my %file2bam;


#smalt map -n 8 -f sam ../raw2/NC_001802_s1k5 A280050.1.fq A280050.2.fq | grep -v \# | 


# 
# 
# 
# Kim Brugger (05 Oct 2012)
sub smalt {

  my $tmp_bam  = EASIH::Pipeline::tmp_file(".bam");
  my $cmd = "smalt map -n 8 -f sam /data/refs/HIV/K03455_s1k6 @$first @$second | egrep -v \\\# | samtools view -t /data/refs/HIV/K03455.fasta.fai -Sb - > $tmp_bam";

  print STDERR "$cmd\n";

  EASIH::Pipeline::submit_job($cmd, "$tmp_bam");

}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_sort {
  my ($input) = @_;

  my $readgroup = $report;
  $readgroup =~ s/\.gz//;
  $readgroup =~ s/\.[1|2].fq//;
  
  my $sample = $readgroup;
  $sample =~ s/\..*//;
  $sample =~ s/_\d*//;

  my $tmp_bam  = EASIH::Pipeline::tmp_file(".bam");
  my $cmd = "$picard -T AddOrReplaceReadGroups.jar I=$input O=$tmp_bam SORT_ORDER=coordinate CN=EASIH PL=$platform LB=$readgroup PU=$run_id  SM=$report VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true";  

  print "$cmd\n";

  EASIH::Pipeline::submit_job($cmd, $tmp_bam);
}




# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub mark_dups {
  my ($input) = @_;

  my $username = scalar getpwuid $<;
  my $metrix_file = EASIH::Pipeline::tmp_file(".mtx");
  my $cmd = "$picard -T MarkDuplicates  I=$input O=$bam_file  M= $metrix_file VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/ MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=true ASSUME_SORTED=Boolean ";
  print STDERR "$cmd\n";
  EASIH::Pipeline::submit_job($cmd, $bam_file);
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
  print "USAGE: $script -1 [fastq file]  -2 [fastq file] -l[oose mapping] -R[eference genome] -d[bsnp rod] -o[ut prefix] -p[latform: illumina or pgm]\n";
  
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
