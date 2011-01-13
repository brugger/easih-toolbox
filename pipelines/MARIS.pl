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

#use lib '/home/cjp64/git/easih-pipeline/modules';
use lib '/home/kb468/easih-pipeline/modules';

use EASIH::JMS;
use EASIH::JMS::Misc;
use EASIH::JMS::Samtools;
use EASIH::JMS::Picard;


our %analysis = ('fastq-split'      => { function   => 'fastq_split',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=02:00:00"},
		 
		 'std-aln'          => { function   => 'bwa_aln',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=12:00:00"},
		 
		 'std-generate'      => { function   => 'bwa_generate',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=3500mb,walltime=12:00:00",},

		 'std-tag_sam'       => { function   => 'sam_add_tags',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},
		 
		 'std-sam2bam'       => { function   => 'EASIH::JMS::Samtools::sam2bam',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=1000mb,walltime=10:00:00"},
		 
		 'std-merge'         => { function   => 'EASIH::JMS::Picard::merge',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=20:00:00",
					  sync       => 1},

		 'std-sort'          => { function   => 'EASIH::JMS::Picard::sort',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=20000mb,walltime=08:00:00"},

		 'std-mark_dup'      =>  { function   => 'EASIH::JMS::Picard::mark_duplicates',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=8000mb,walltime=16:00:00"},

		 'std-index'         => { function   => 'EASIH::JMS::Samtools::index',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
		 

		 'get_all_mapped'    => { function   => 'EASIH::JMS::Samtools::get_mapped',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'get_all_unmapped'  => { function   => 'EASIH::JMS::Samtools::get_unmapped',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'identify_indel'    => { function   => 'identify_indel',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'realign_indel'     => { function   => 'realign_indel',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=08:00:00"},

		 'realigned_merge'   => { function   => 'EASIH::JMS::Picard::merge',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=20:00:00",
					  sync       => 1},

		 'realigned_rename'  => { function   => 'rename' },

		 
		 'realigned_sort'    => { function   => 'EASIH::JMS::Picard::sort',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=20000mb,walltime=08:00:00"},
		 

		 'realigned_index'   => { function   => 'EASIH::JMS::Samtools::index',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
		 
		 'call_indels'       => { function   => 'call_indels',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},
		 
		 'merge_indels'      => { function   => 'merge_indels',
					  sync       => 1},
		 

		 'identify_snps'     => { function   => 'identify_snps',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'filter_snps'       => { function   => 'filter_snps',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500b,walltime=01:00:00"},
		 

		 'merge_vcfs'        => { function   => 'merge_vcfs',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=02:00:00", 
					  sync       => 1},
		 
		 'cluster_snps'      => { function   => 'cluster_snps',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=01:00:00"},
		 
		 'rescore_snps'      => { function   => 'rescore_snps',
					  hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=50000mb,walltime=01:00:00"},		 
		 
		 

    );
		     

our %flow = ( 'csfasta2fastq'     => 'std-aln',
	      'fastq-split'       => 'std-aln',

	      'std-aln'           => 'std-generate',
	      'std-generate'      => 'std-tag_sam',
	      'std-tag_sam'       => 'std-sam2bam',
	      'std-sam2bam'       => 'std-merge',
	      'std-merge'         => 'std-sort',
	      'std-sort'          => 'std-index',
	      'std-mark_dup'      => 'std-index',

	      'std-index'        => ['identify_indel', 'get_all_unmapped'],
	      'get_all_unmapped' => 'realigned_merge',
	      'identify_indel'   => 'realign_indel',
	      'realign_indel'    => 'realigned_merge',
	      'realigned_merge'  => 'realigned_sort',
	      'realigned_sort'   => 'realigned_rename',
	      'realigned_rename' => 'realigned_index',

	      'realigned_index'  => ['call_indels','identify_snps'],
	      'call_indels'      => 'merge_indels',

	      'identify_snps'    => 'filter_snps',
	      'filter_snps'      => 'merge_vcfs',
	      'merge_vcfs'       => 'cluster_snps',
	      'cluster_snps'     => 'rescore_snps'
	      );

#EASIH::JMS::no_store();
#EASIH::JMS::print_flow('fastq-split');

my %opts;
getopts('1:2:d:e:f:hH:I:lmM:n:No:p:Q:Pr:R:sS:L:v', \%opts);

usage() if ( $opts{h});
my $hard_reset    = $opts{'H'};
my $soft_reset    = $opts{'S'};

if ( $soft_reset ) {
  print "Doing a soft reset/restart\n";
  &EASIH::JMS::reset($soft_reset);
  getopts('1:2:d:e:f:hH:lmM:n:No:p:Q:Pr:R:S:L:v', \%opts);
}
elsif ( $hard_reset ) {
  &EASIH::JMS::hard_reset($hard_reset);
  getopts('1:2:d:e:f:hH:lmM:n:No:p:Q:Pr:R:S:L:v', \%opts);
}



my $username = scalar getpwuid $<;

# if using standard naming, this is a lot easier.
if ( $opts{Q} ) {
  $opts{'1'} = "$opts{Q}.1.fq"    if ( -e "$opts{Q}.1.fq");
  $opts{'1'} = "$opts{Q}.1.fq.gz" if ( -e "$opts{Q}.1.fq.gz");
  $opts{'2'} = "$opts{Q}.2.fq"    if ( -e "$opts{Q}.2.fq");
  $opts{'2'} = "$opts{Q}.2.fq.gz" if ( -e "$opts{Q}.2.fq.gz");
  $opts{'L'} = "$opts{Q}.log";
  $opts{'o'} = "$opts{Q}";
  $opts{'l'} = 1;
}  

my $first         = $opts{'1'}     || usage();
my $second        = $opts{'2'};
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
my $report        = $opts{'o'}     || usage();
my $platform      = uc($opts{'p'}) || usage();
$platform = 'SOLEXA'      if ( $platform eq 'ILLUMINA');
my $print_filter  = $opts{'P'};
my $readgroup     = $opts{'r'} || $report;
my $reference     = $opts{'R'}     || usage();
my $no_sw_pair    = $opts{'s'};
my $align_param   = ' ';

my $bam_file      = "$report.bam";

#EASIH::JMS::verbosity(100) if ( $opts{v});



open (*STDOUT, ">> $log") || die "Could not open '$log': $!\n" if ( $log );


# set platform specific bwa aln parameters
$align_param .= " -c "      if ( $platform eq "SOLID");
$align_param .= " -q 15 "   if ( $platform eq "SOLEXA");
# and loose mapping
$align_param .= " -e5 "     if ( $loose_mapping);

$no_sw_pair     = 1 if ($platform eq "SOLID");
my $sampe_param = "";
$sampe_param    = '-s ' if ( $first && $second && $no_sw_pair);
$sampe_param    = "-M $insert_size "  if ( $first && $second && $insert_size);




# Only paired ends runs gets marked duplicates.
$flow{'std-sort'} = 'std-mark_dup' if (($first && $second) || $mark_dup );

my $bwa          = EASIH::JMS::Misc::find_program('bwa');
my $fq_split     = EASIH::JMS::Misc::find_program('fastq_split.pl');
my $samtools     = EASIH::JMS::Misc::find_program('samtools');
my $tag_sam      = EASIH::JMS::Misc::find_program('tag_sam.pl');
my $gatk         = EASIH::JMS::Misc::find_program('gatk');

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
$extra_report .= "easih-pipeline: " . EASIH::JMS::version() . "\n";

$extra_report .= "align_param ==> $align_param and $sampe_param\n";
$extra_report .= "Binaries used..\n";
$extra_report .= "BWA: " . EASIH::JMS::Misc::bwa_version( $bwa ) . "\n";
$extra_report .= "Samtools: " . EASIH::JMS::Samtools->version() ."\n";
$extra_report .= "GATK: " .`$gatk --version`;
$extra_report .= "Picard: " . EASIH::JMS::Picard->version() ."\n";
$extra_report .= "Command line: $0 ".EASIH::JMS::args() ."\n";


EASIH::JMS::mail_report($email, $bam_file, $extra_report);


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
    $cmd = "$bwa sampe $sampe_param  $reference $$input{first_sai} $$input{second_sai} $$input{first_fq} $$input{second_fq} > $tmp_file";
  }
  else {
    $cmd = "$bwa samse -f $tmp_file $reference $$input{first_sai} $$input{first_fq}";
  }

  EASIH::JMS::submit_job($cmd, $tmp_file);
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
    my $cmd = "$gatk -T RealignerTargetCreator -R $reference -o $tmp_file -I $input -L $name";
    EASIH::JMS::submit_job($cmd, "$tmp_file $name $input");
  }
  
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
    $cmd = "$gatk -T IndelRealigner -targetIntervals $interval_file -L $region --output $tmp_file -R $reference -I $tmp_bam_file";
  }

  EASIH::JMS::submit_job($cmd, $tmp_file);
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

  my $cmd = "cat  @inputs > $report";
  EASIH::JMS::submit_system_job("cat  @inputs > $report.indels");
}


sub rename {
  my ($input) = @_;

  EASIH::JMS::submit_system_job("mv $input $bam_file", $bam_file);
}


sub identify_snps {

  my @names = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }

  foreach my $name ( @names ) {
    my $tmp_file = EASIH::JMS::tmp_file(".vcf");
    my $cmd = "$gatk -T UnifiedGenotyper -R $reference -I $bam_file -G Standard -D $dbsnp -varout $tmp_file -L $name ";
    $cmd .= " -pl $platform " if ( $platform);
    EASIH::JMS::submit_job($cmd, "$tmp_file  -L $name");
  }
  
}

sub filter_snps {
  my ($input) = @_;

  my ($input_file, $region) = split(" ", $input);

  my $tmp_file = EASIH::JMS::tmp_file(".filtered.vcf");

  my $entries = `egrep -cv \# $input_file`;
  chomp( $entries );
  # Filtering will fail on an empty file, so we just fake it for now
  # and ensure that the pipeline will not break its dependencies for a restart.
  if ( $entries == 0 ) {
    EASIH::JMS::submit_system_job("touch $tmp_file", $tmp_file);
  }
  else {
    my $cmd = "$gatk -T VariantFiltration  -R $reference  -B variant,VCF,$input  -o $tmp_file $filter";
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}		



sub merge_vcfs {
  my (@inputs) = @_;

  my $merged_file = EASIH::JMS::tmp_file(".merged.vcf");
  
  if ( @inputs == 1 ) {
    EASIH::JMS::submit_system_job("mv @inputs $merged_file", $merged_file);
  }
  else {
    my $cmd = "$gatk -T CombineVariants -R $reference -o $merged_file -variantMergeOptions UNION -genotypeMergeOptions UNIQUIFY";
    my $count = 1;
    foreach my $input ( @inputs ) {
      $cmd .= " -B variant,VCF,$input ";
      $count++;
    }
    EASIH::JMS::submit_job($cmd, $merged_file);
  }
}

sub cluster_snps {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".cluster");
  my $cmd = "$gatk -T GenerateVariantClusters -R $reference -B input,VCF,$input --DBSNP $dbsnp -an QD -an SB -an HaplotypeScore -an HRun -clusterFile $tmp_file";
  EASIH::JMS::submit_system_job($cmd, "$tmp_file -B input,VCF,$input");
}		
		
		

sub rescore_snps {
  my ($input) = @_;

  $report =~ s/.vcf\z//;
  my $cmd = "$gatk -T VariantRecalibrator -R $reference --DBSNP $dbsnp -clusterFile $input -output  $report.snps --target_titv 3.0 ";
  EASIH::JMS::submit_system_job($cmd, $report);
}		




# 
# Ensure that the reference and it auxiliary files are all present.
# 
# Kim Brugger (02 Aug 2010)
sub validate_input {
  
  my @errors;
  my @info;

  push @errors, "$first does not exist"  if ( ! -e $first );
  push @errors, "$second does not exist"  if ( $second &&  ! -e $second );


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
  push @errors, "'$dbsnp' does end with .rod as expected\n" if ($dbsnp !~ /.rod\z/);

  push @errors, "Platform must be either SOLEXA or SOLID not '$platform'" if ( $platform ne "SOLEXA" && $platform ne 'SOLID');

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
  print "extra flags: -s[ disable Smith-Waterman for the unmapped mate]\n";
  print "extra flags: -S[oft reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "\n";

  print "easih-pipeline: " . &EASIH::JMS::version() . "\n";

  exit;

}