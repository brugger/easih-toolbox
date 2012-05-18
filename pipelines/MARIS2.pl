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

use lib '/home/kb468/easih-pipeline/modules';
use lib '/home/kb468/easih-toolbox/modules';

use EASIH::JMS;
use EASIH::JMS::Misc;
use EASIH::JMS::Samtools;
use EASIH::JMS::Picard;


our %analysis = ('fastq-split'      => { function   => 'fastq_split',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=02:05:00"},
		 
		 'std-aln'          => { function   => 'bwa_aln',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=4,walltime=08:05:00"},
		 
		 'sai2bam'          => { function   => 'bwa_sai2bam',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=4,walltime=04:15:00",},

		 'std-sort'          => { function   => 'EASIH::JMS::Picard::sort',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=05:05:00"},

		 'std-merge'         => { function   => 'EASIH::JMS::Picard::merge',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=08:00:00",
					  sync       => 1},

		 'std-mark_dup'      => { function   => 'EASIH::JMS::Picard::mark_duplicates',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=16:00:00"},
		 
		 'std-index'         => { function   => 'EASIH::JMS::Samtools::index',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=04:00:00"},
		 
		 'get_all_mapped'    => { function   => 'EASIH::JMS::Samtools::get_mapped',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=08:00:00"},
		 
		 'get_all_unmapped'  => { function   => 'EASIH::JMS::Samtools::get_unmapped',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=08:00:00"},
		 
		 'identify_indel'    => { function   => 'identify_indel_mpi',
					  hpc_param  => "-NEP-fqs -l nodes=2:ppn=4,walltime=02:00:00"},
		 
		 'realign_indel'     => { function   => 'realign_indel_mpi',
					  hpc_param  => "-NEP-fqs -l nodes=2:ppn=4,walltime=08:00:00"},

		 'realigned_merge'   => { function   => 'EASIH::JMS::Picard::merge',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=20:00:00",
					  sync       => 1},

		 'realigned_sort'    => { function   => 'EASIH::JMS::Picard::sort',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=08:00:00"},
		 
		 'realigned_rename'  => { function   => 'rename' },

		 'realigned_index'   => { function   => 'EASIH::JMS::Samtools::index',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=04:00:00"},
		 
		 'identify_variation'     => { function   => 'identify_variation_mpi',
					  hpc_param  => "-NEP-fqs -l nodes=3:ppn=4,walltime=02:00:00"},
		 
		 'filter_variation'       => { function   => 'filter_variation',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=02:00:00"},
		 

		 'merge_vcfs'        => { function   => 'merge_vcfs',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=02:00:00"},
		 
		 'cluster_snps'      => { function   => 'cluster_snps',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=01:00:00"},
		 
		 'rescore_snps'      => { function   => 'rescore_snps',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,walltime=01:00:00"},

		 'stats'             => { function   => 'stats',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=4,walltime=02:00:00"},

		 'finished'          => { function   => 'finished',
					  },
    );
		     

our %flow = ( 
	      'fastq-split'       => 'std-aln',

	      'std-aln'           => 'sai2bam',
	      'sai2bam'           => 'std-sort',
	      'std-sort'          => 'std-merge',
              'std-merge'         => 'std-index',
              'std-mark_dup'      => 'std-index',

	      'std-index'        => ['identify_indel', 'get_all_unmapped'],
	      'get_all_unmapped' => 'realigned_merge',
	      'identify_indel'   => 'realign_indel',
	      'realign_indel'    => 'realigned_merge',
	      'realigned_merge'  => 'realigned_rename',
	      'realigned_rename' => 'realigned_index',
	      
	      'realigned_index'  => 'identify_variation',

	      'identify_variation'=> 'merge_vcfs',
	      'merge_vcfs'       =>  'filter_variation',
#	      'merge_vcfs'       =>  ['get_indels','rescore_snps'],
#	      'get_indels'       => 'filter_variation',
#	      'rescore_snps'     => 'filter_variation',
	      'filter_variation' => 'stats',
	      'stats'            => 'finished'
	      );

#EASIH::JMS::no_store();
#EASIH::JMS::print_flow('fastq-split');

my $opts = '1:2:d:De:f:hH:I:lL:mM:n:No:p:Q:Pr:R:sS:vV';
my %opts;
getopts($opts, \%opts);

usage() if ( $opts{h});
my $hard_reset    = $opts{'H'};
my $soft_reset    = $opts{'S'};

if ( $soft_reset ) {
  print "Doing a soft reset/restart\n";
  &EASIH::JMS::reset($soft_reset);
  getopts($opts, \%opts);
}
elsif ( $hard_reset ) {
  print "Doing a hard reset/restart\n";
  &EASIH::JMS::hard_reset($hard_reset);
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

  if ( $opts{'r'} ) {
    $opts{'r'} = 1;
  }

  my $ref_dir       = $opts{'r'};


  my $freeze_file = "$opts{'o'}.maris";
  system "mv $freeze_file $freeze_file.backup"  if ( -e $freeze_file );

  EASIH::JMS::freeze_file($freeze_file);
}  

my $first         = $opts{'1'}     || usage();
$first            = [split(",", $first)];
my $second        = $opts{'2'};
$second           = [split(",", $second)] if ( $second );
my $dbsnp         = $opts{'d'}     || usage();
my $email         = $opts{'e'}     || "$username\@cam.ac.uk";
my $filter        = $opts{'f'}     || "exon";
my $insert_size   = $opts{'I'};
my $loose_mapping = $opts{'l'}     || 0;
my $log           = $opts{'L'};
my $mark_dup      = $opts{'m'};
my $min_depth     = $opts{'M'}     || 0;
my $split         = $opts{'n'}     || 5000000;
my $no_split      = $opts{'N'}     || 0;
our $report       = $opts{'o'}     || usage();
my $platform      = uc($opts{'p'}) || usage();
$platform = 'ILLUMINA'      if ( $platform eq 'ILLUMINA');
my $print_filter  = $opts{'P'};
our $reference    = $opts{'R'}     || usage();
my $no_sw_pair    = $opts{'D'};
my $align_param   = ' ';

my $bam_file      = "$report.bam";

#EASIH::JMS::verbosity(100) if ( $opts{v});


my $scrub_data    = $opts{s} || 0;
# change the flow to accommodate scrubbing...
if ( $scrub_data && $platform eq "SOLID" ) {

  $flow{'realigned_sort'}   = 'scrub_index';
  $flow{'scrub_index'}      = 'scrub';
  $flow{'scrub'}            = 'scrub_rename';
  $flow{'scrub_rename'}     = 'realigned_index';
  push @{$flow{'realigned_index'}}, 'scrub_bed';
}


open (*STDOUT, ">> $log") || die "Could not open '$log': $!\n" if ( $log );


# set platform specific bwa aln parameters
$align_param .= " -c "      if ( $platform eq "SOLID");
$align_param .= " -q 15 "   if ( $platform eq "ILLUMINA");
# and loose mapping
$align_param .= " -e5 "     if ( $loose_mapping);

$no_sw_pair     = 1 if ($platform eq "SOLID");
my $sampe_param = " -P ";
$sampe_param    .= '-s ' if ( $first && $second && $no_sw_pair);
$sampe_param    .= "-M $insert_size "  if ( $first && $second && $insert_size);


# Only paired ends runs gets marked duplicates.
$flow{'std-merge'} = 'std-mark_dup' if (($first && $second) || $mark_dup );

#$flow{'std-merge'} = 'std-index' if (($first && $second) || $mark_dup );


my $bwa             = EASIH::JMS::Misc::find_program('bwa');
my $fq_split        = EASIH::JMS::Misc::find_program('fastq_split.pl');
my $samtools        = EASIH::JMS::Misc::find_program('samtools');
my $tag_sam         = EASIH::JMS::Misc::find_program('tag_sam.pl');
my $gatk            = EASIH::JMS::Misc::find_program('gatk');
my $scrubber        = EASIH::JMS::Misc::find_program('SOLiD_scrubber.pl');
my $bam2scubber_bed = EASIH::JMS::Misc::find_program('bam2scubber-bed.pl');

$samtools = "/home/kb468/bin/samtools";
$gatk     = "/home/kb468/bin/gatk112";
$gatk     = "/home/kb468/bin/gatk1314";
my $mpi_q = "/home/kb468/easih-pipeline/tools/mpiexec_queue.pl";

validate_input();


if ($filter eq "wgs" ) {
  $min_depth ||= 20;

  $filter  = "--filterExpression 'DP < $min_depth' --filterName shallow$min_depth ";
#  $filter .= "-clusterWindowSize 10 ";
  $filter .= "--filterExpression 'AB > 0.75 && DP > 40 || DP > 100 || MQ0 > 40 || SB > -0.10'   --filterName StdFilter ";
  $filter .= "--filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)'  --filterName HARD_TO_VALIDATE ";
  
}
elsif ( $filter eq "wgs-low" ) {
  $min_depth ||= 5;

  $filter  = "--filterExpression 'DP < $min_depth' --filterName shallow$min_depth ";
  $filter .= "--filterExpression 'MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1'  --filterName HARD_TO_VALIDATE";
}
elsif ( $filter eq "exon" ) {
  $min_depth ||= 20;

  $filter  = "--filterExpression 'DP < $min_depth' --filterName shallow$min_depth ";
  $filter .= " --filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' -filterName StddFilter ";
  $filter .= " --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE";  
}
elsif ( $filter eq "exon-low" ) {
  $min_depth ||= 5;

  $filter  = "--filterExpression 'DP < $min_depth' --filterName shallow$min_depth ";
  $filter .= " --filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' -filterName StdFilter ";
  $filter .= " --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE";  
}

if ( $print_filter ) {
  print "GATK Filter to be used: $filter\n";
  exit;
}

#EASIH::JMS::verbosity(10);
EASIH::JMS::backend('Darwin');
#EASIH::JMS::backend('Local');
EASIH::JMS::max_retry(0);


if ( $opts{ V } ) {

  &EASIH::JMS::run('identify_variation');
  
  
}
elsif ( $no_split ) {
  &EASIH::JMS::run('std-aln');
}
else {
  &EASIH::JMS::run('fastq-split');
}

&EASIH::JMS::store_state();

my $extra_report = "1 ==> @$first\n";
$extra_report .= "2 ==> @$second\n" if ( @$second );
$extra_report .= "bamfile ==> $bam_file\n";
$extra_report .= "snp_file ==> $report.snps\n";
$extra_report .= "indel_file ==> $report.indel\n";
$extra_report .= "easih-pipeline: " . EASIH::JMS::version() . "\n";

$extra_report .= "align_param ==> $align_param and $sampe_param\n";
$extra_report .= "Binaries used..\n";
$extra_report .= "BWA: " . EASIH::JMS::Misc::bwa_version( $bwa ) . "\n";
$extra_report .= "Samtools: " . EASIH::JMS::Samtools->version() ."\n";
#$extra_report .= "GATK: " .`$gatk --version`;
$extra_report .= "Picard: " . EASIH::JMS::Picard->version() ."\n";
$extra_report .= "Command line: $0 ".EASIH::JMS::args() ."\n";


EASIH::JMS::mail_report($email, $bam_file, $extra_report);

#EASIH::JMS::delete_tmp_files();

my %split2files;

sub fastq_split {

  foreach my $first_file ( @$first ) {
    my $second_file = shift @$second if ( @$second );

    my $tmp_file = EASIH::JMS::tmp_file();
    # keep track of the original file
    $split2files{$tmp_file} = $first_file; 
    
    my $cmd = "$fq_split -e $split  -1 $first_file";
    $cmd .= " -2 $second_file "  if ( $second_file);
    $cmd .= " -o tmp/ > $tmp_file";
    
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}



sub bwa_aln {
  my ($input) = @_;

  if ( $no_split ) {
    
    if ( $first && $second ) {
      foreach my $first_file ( @$first ) {
	my $second_file = shift @$second if ( @$second );
	$split2files{$first_file} = $first_file; 

	my $first_tmp_file  = EASIH::JMS::tmp_file(".sai");
	my $second_tmp_file = EASIH::JMS::tmp_file(".sai");
	my $cmd = "$bwa aln -t 4 $align_param  -f $first_tmp_file  $reference $first_file ;";
	$cmd   .= "$bwa aln -t 4 $align_param  -f $second_tmp_file $reference $second_file ";
	
	my $output = { "first_fq"   => $first_file,
		       "first_sai"  => $first_tmp_file,
		       "second_fq"  => $second_file,
		       "second_sai" => $second_tmp_file};
	
#      EASIH::JMS::submit_job($cmd, "$first_tmp_file $second_tmp_file $first $second");
	EASIH::JMS::submit_job($cmd, $output);
      }
    }
    else {
      foreach my $first_file ( @$first ) {
	$split2files{$first_file} = $first_file; 
	my $tmp_file  = EASIH::JMS::tmp_file(".sai");
	my $cmd = "$bwa aln -t 4 $align_param  -f $tmp_file $reference $first_file ";
	my $output = { "first_fq"   => $first_file,
		       "first_sai"  => $tmp_file};
	EASIH::JMS::submit_job($cmd, $output);
      }
    }
  }
  else {

    open (my $files, $input) || die "Could not open '$input': $!\n";
    while (<$files>) {
      chomp;
      my ($file1, $file2) = split("\t", $_);

      $split2files{ $file1 } = $split2files{ $input };
      $split2files{ $file2 } = $split2files{ $input } if ( $file2 );
      
      if ( $file1 && $file2 ) {
	my $first_tmp_file  = EASIH::JMS::tmp_file(".sai");
	my $second_tmp_file = EASIH::JMS::tmp_file(".sai");
	my $cmd = "$bwa aln -t 4 $align_param  -f $first_tmp_file  $reference $file1 ;";
	$cmd   .= "$bwa aln -t 4 $align_param  -f $second_tmp_file $reference $file2 ";
	my $output = { "first_fq"   => $file1,
		       "first_sai"  => $first_tmp_file,
		       "second_fq"  => $file2,
		       "second_sai" => $second_tmp_file};
		     
	EASIH::JMS::submit_job($cmd, $output);
      }
      else {
	my $tmp_file  = EASIH::JMS::tmp_file(".sai");
	my $cmd = "$bwa aln -t 4 $align_param  -f $tmp_file $reference $file1 ";
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
    $cmd = "$bwa sampe $sampe_param  $reference $$input{first_sai} $$input{second_sai} $$input{first_fq} $$input{second_fq} > $tmp_file";
  }
  else {
    $cmd = "$bwa samse -f $tmp_file $reference $$input{first_sai} $$input{first_fq}";
  }

  EASIH::JMS::submit_job($cmd, $tmp_file);
}



# 
# 
# 
# Kim Brugger (01 Aug 2011)
sub bwa_sai2bam {

  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".bam");

  my $cmd;
  if (defined($$input{'second_sai'}) ) {
    $cmd = "$bwa sampe $sampe_param  $reference $$input{first_sai} $$input{second_sai} $$input{first_fq} $$input{second_fq} ";
  }
  else {
    $cmd = "$bwa samse  $reference $$input{first_sai} $$input{first_fq} ";
  }

  my $readgroup = $split2files{ $$input{first_fq} };
  $readgroup =~ s/\.gz//;
  $readgroup =~ s/\.[1|2].fq//;

  
  my $sample = $readgroup;
  $sample =~ s/\..*//;
  $sample =~ s/_\d*//;

  
  $cmd .= " | $tag_sam -R $reference -p $platform -r $readgroup -s $sample ";
  $cmd .= " | $samtools view -T $reference -Sb -  >  $tmp_file ";

  EASIH::JMS::submit_job($cmd, $tmp_file);

  
}


sub sortNcalmd { 
  my ($input) = @_;

  my $username = scalar getpwuid $<;
  my $picard = EASIH::JMS::Misc::find_program('picard');

  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd = "$picard -T SortSam  I= $input O= /dev/stdout SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/ | samtools calmd -eb - $reference > $tmp_file";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}



sub calmd {
  my ($input) = @_;
  
  my $tmp_file = EASIH::JMS::tmp_file(".calmd.bam");
  my $cmd = "$samtools calmd -b $input $reference > $tmp_file";


  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub calmd_mpi {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".calmd.bam");
  my $cmd = "$samtools calmd -b $input $reference > $tmp_file";


  EASIH::JMS::submit_job($cmd, $tmp_file);
}




sub stats {
  my ($input) = @_;
  
  my $commands_file = EASIH::JMS::tmp_file(".mpiexec");
  open(my $cmds, ">$commands_file") || die "Could not write to '$commands_file': $!\n";

  print  $cmds "$samtools flagstat $report.bam > $report.flagstat\n";
  print  $cmds "md5sum $report.bam > $report.bam.md5 \n";
  print  $cmds "md5sum $report.bam.bai > $report.bam.bai.md5 \n";
  print  $cmds "md5sum $report.vcf > $report.vcf.md5 \n";
  print  $cmds "md5sum $report.maris > $report.maris.md5\n";

  close ($cmds );

  my $cmd = "$mpi_q -c 4 $commands_file ";
  EASIH::JMS::submit_job($cmd, "");
}



sub identify_indel {
  my ( $input ) = @_;

  my @names = ();
  open(my $spipe, "$samtools view -H $input | ") || die "Could not open '$input': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }

  foreach my $name ( @names ) {
    my $tmp_file = EASIH::JMS::tmp_file(".intervals");
    my $cmd = "$gatk -T RealignerTargetCreator -R $reference -o $tmp_file  -L $name -I $input";
    EASIH::JMS::submit_job($cmd, "$tmp_file $name $input");
  }
  
}


sub identify_indel_mpi {
  my ( $input ) = @_;

  my @names = ();
  open(my $spipe, "$samtools view -H $input | ") || die "Could not open '$input': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }


  my @targets;
  
  my $commands_file = EASIH::JMS::tmp_file(".mpiexec");
  open(my $cmds, ">$commands_file") || die "Could not write to '$commands_file': $!\n";

  foreach my $name ( @names ) {
    my $tmp_file = EASIH::JMS::tmp_file(".intervals");
    push @targets, "$tmp_file  -L $name -I $input";
    print $cmds "$gatk -T RealignerTargetCreator -R $reference -o $tmp_file  -L $name -I $input\n";
  }
  close ($cmds );

  my $cmd = "$mpi_q -c 8 $commands_file ";
  EASIH::JMS::submit_job($cmd, \@targets);
  
}


sub realign_indel_mpi {
  my ($inputs) = @_;

  my @targets;
  my $commands_file = EASIH::JMS::tmp_file(".mpiexec");
  open(my $cmds, ">$commands_file") || die "Could not write to '$commands_file': $!\n";

  foreach my $input ( @$inputs ) {

    my $tmp_file = EASIH::JMS::tmp_file(".bam");


  # If the interval file is empty the realigner ignores the region and produces an empty bamfile...
    if ( 0 ) {# -z $interval_file ) {
#      print $cmds "$samtools view -b $tmp_bam_file $region > $tmp_file";
    }
    else { #
      print $cmds "$gatk -T IndelRealigner  -targetIntervals $input -o $tmp_file -R $reference -baq CALCULATE_AS_NECESSARY\n";
      push @targets, $tmp_file;
    }
  }

  close($cmds);

  my $cmd = "$mpi_q -c 8 $commands_file ";
  EASIH::JMS::submit_job($cmd, \@targets);
}

sub realign_indel {
  my ($input) = @_;

  my ($interval_file, $region, $tmp_bam_file) = split(" ", $input);

  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd;
  # If the interval file is empty the realigner ignores the region and produces an empty bamfile...
  if (  -z $interval_file ) {
    $cmd = "$samtools view -b $tmp_bam_file $region > $tmp_file";
  }
  else {
    $cmd = "$gatk -T IndelRealigner -baq CALCULATE_AS_NECESSARY -targetIntervals $interval_file -L $region -o $tmp_file -R $reference -I $tmp_bam_file ";
  }

  EASIH::JMS::submit_job($cmd, $tmp_file);
}



# 
# 
# 
# Kim Brugger (18 Jan 2011)
sub solid_scrubber {
  my ($input) = @_;
  
  my $tmp_file = EASIH::JMS::tmp_file(".scrubbed.bam");
  my $cmd = "$scrubber -b $input -R $reference -Ul $report.scrubber.txt  > $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);

}


# 
# 
# 
# Kim Brugger (18 Jan 2011)
sub scrub_bed {
  my ($input) = @_;

  my $outfile = "$report\_scrub.bed";

  my $cmd = "$bam2scubber_bed  $input > $outfile";

  EASIH::JMS::submit_job($cmd, $outfile);

}



# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub merge_indels {
  my (@inputs) = @_;

  my $cmd = "cat  @inputs > $report";
  EASIH::JMS::submit_system_job("cat  @inputs > $report.indels.vcf");
}


sub rename {
  my ($input) = @_;

  EASIH::JMS::submit_system_job("mv $input $bam_file", $bam_file);
}


sub identify_variation_mpi {

  my @names = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }


  my @targets;

  my $commands_file = EASIH::JMS::tmp_file(".mpiexec");
  open(my $cmds, ">$commands_file") || die "Could not write to '$commands_file': $!\n";

  foreach my $name ( @names ) {
    my $tmp_file = EASIH::JMS::tmp_file(".vcf");
    my $cmd = "$gatk -T UnifiedGenotyper -R $reference -I $bam_file -G Standard -D $dbsnp -o $tmp_file -L $name -glm BOTH ";
    print $cmds "$cmd\n";
    push @targets, "$tmp_file";

  }
  close ($cmds );

  my $cmd = "$mpi_q -c 12 $commands_file ";
  EASIH::JMS::submit_job($cmd, \@targets);
}

sub filter_variation {
  my ($input) = @_;

  my $merged_file = "$report.vcf";

  my $entries = `egrep -cv \# $input`;
  chomp( $entries );
  # Filtering will fail on an empty file, so we just fake it for now
  # and ensure that the pipeline will not break its dependencies on a restart.
  if ( $entries == 0 ) {
    EASIH::JMS::submit_system_job("mv $input $merged_file", "");
  }
  else {
    my $cmd = "$gatk -T VariantFiltration  -R $reference  -V $input  -o $merged_file $filter ";
    EASIH::JMS::submit_job($cmd, "$merged_file");
  }

}

sub merge_vcfs {
  my ($inputs) = @_;

#  my $merged_file = EASIH::JMS::tmp_file(".merged.vcf");
  my $tmp_file = EASIH::JMS::tmp_file(".merged.vcf");
  
  if ( @$inputs == 1 ) {
    EASIH::JMS::submit_system_job("mv @$inputs $tmp_file", $tmp_file);
  }
  else {
    my $cmd = "$gatk -T CombineVariants -R $reference -o $tmp_file ";
    my $count = 1;
    foreach my $input ( @$inputs ) {
      $cmd .= " -V:variant $input ";
      $count++;
    }
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}


# 
# Pulls out all the indels, so we can merge them back in later.
# 
# Kim Brugger (02 Dec 2011)
sub get_indels {
  my ($input) = @_;

  my $tmp_indel_file = EASIH::JMS::tmp_file(".indels.vcf");
  my $cmd = "$gatk -T SelectVariants  -indels -R $reference -V $input -o $tmp_indel_file";

  

}



# 
# This only works for human data...
# 
# Kim Brugger (02 Dec 2011)
sub rescore_snps {
  my ($input) = @_;


  my $tmp_cmd_file   = EASIH::JMS::tmp_file(".sh");
  open(my $cmds, "> $tmp_cmd_file") || die "Could not write to '$tmp_cmd_file': $!\n";

  my $tmp_snp_file      = EASIH::JMS::tmp_file(".snps.vcf");
  my $tmp_recalc_file   = EASIH::JMS::tmp_file(".recalc");
  my $tmp_tranches_file = EASIH::JMS::tmp_file(".tranches");
  my $tmp_recal_snps    = EASIH::JMS::tmp_file(".rc.snps.vcf");

  # get all the snps from the VCF file...
  print  $cmds "$gatk -T SelectVariants  -snps   -R $reference -V $input -o $tmp_snp_file\n";
  #calculate the recalibration thingy


#ref='/home/easih/refs/human_1kg/human_g1k_v37.fasta'
#dbsnp='/home/easih/refs/GATK/dbsnp_132.b37.vcf'
#indels='/home/easih/refs/GATK/1000G_indels_for_realignment.b37.vcf'
#omni='/home/easih/refs/GATK/1000G_omni2.5.b37.sites.vcf'
#hapmap='/home/easih/refs/GATK/hapmap_3.3.b37.sites.vcf'


#  print  $cmds "$gatk -T VariantRecalibrator -nt 4  -R $reference -input $tmp_snp_file ";
#  print  $cmds "-resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $hapmap ";
#  print  $cmds "-resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 $omni "; 
#  print  $cmds "-resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $dbsnp "; 
#  print  $cmds "-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode SNP -mG 6 -recalFile $tmp_recalc_file -tranchesFile $tmp_tranches_file -rscriptFile ${name}_var_recal.plots.R\n";

 print $cmds "$gatk -T ApplyRecalibration -R $reference -B:input,VCF $tmp_snp_file -recalFile $tmp_recalc_file -tranchesFile $tmp_tranches_file -o $tmp_recal_snps\n";


  
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
  print $out "$report.vcf\n";
  print $out "$report.vcf.md5\n";
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


  push @errors, "'$dbsnp' does not exists\n" if (! -e $dbsnp);
  push @errors, "'$dbsnp' don't end with .vcf as expected\n" if ($dbsnp !~ /.vcf\z/);

  push @errors, "Platform must be either ILLUMINA or SOLID not '$platform'" if ( $platform ne "ILLUMINA" && $platform ne 'SOLID');

  push @errors, "$filter should be one of the following: wgs,wgs-low,exon,exon-low\n"
      if ($filter ne "wgs" && $filter ne "wgs-low" && $filter ne "exon" && $filter ne "exon-low");

  push @errors, "-M '$min_depth' is <0 or not a number\n" if ($min_depth < 0 || $min_depth !~ /^\d+$/);

  # print the messages and die if critical ones.
  die join("\n", @errors) . "\n"   if ( @errors );
  print  join("\n", @info) . "\n"   if ( @info );
}



# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {

  my $script = $0;
  $script =~ s/.*\///;
  print "USAGE: $script -1 [fastq file]  -2 [fastq file] -l[oose mapping] -R[eference genome] -d[bsnp rod] -o[ut prefix] -p[latform: illumina or solid]\n";
  
  print "\nor extrapolate the standard <fq, log, out names> with the -Q flag\n";
  print "EXAMPLE: $script -Q [base name] -l[oose mapping] -R[eference genome] -d[bsnp rod] -p[latform: illumina or solid]\n";

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
  print "extra flags: -P[rint GATK filter, and exit]\n"; 
  print "extra flags: -r[ead group]\n"; 
  print "extra flags: -s[crub mapping (Solid only!)]\n";
  print "extra flags: -S[oft reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "\n";

  print "easih-pipeline: " . &EASIH::JMS::version() . "\n";

  use EASIH::Toolbox;
  print "easih-toolbox: " . &EASIH::Toolbox::version() . "\n";

  exit;

}
