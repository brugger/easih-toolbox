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



our %analysis = ('csfasta2fastq'    => { function   => 'csfasta2fastq',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=05:00:00"},

		 'fastq-split'      => { function   => 'fastq_split',
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
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00",
					 sync       => 1},

		 'get_mapped'        => { function   => 'EASIH::JMS::Samtools::get_mapped',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'get_unmapped'      => { function   => 'EASIH::JMS::Samtools::get_unmapped',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'unmapped_bam2fq'   => { function   => 'bam2fq',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=01:00:00"},

		 '2nd-aln'           => { function   => 'bwa_aln_loose',
				       hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=12:00:00"},
		 
		 '2nd-generate'      => { function   => 'bwa_generate',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=10:00:00",},

		 '2nd-sam2bam'       => { function   => 'EASIH::JMS::Samtools::sam2bam',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=1000mb,walltime=10:00:00"},
		 
		 'mapped_merge'      => { function   => 'EASIH::JMS::Picard::merge',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00",
					 sync       => 1},

		 'mapped_sort'       => { function   => 'EASIH::JMS::Picard::sort',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=20000mb,walltime=08:00:00"},


		 'mapped_index'    =>{ function   => 'EASIH::JMS::Samtools::index',
					hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
		 


		 'get_all_mapped'   => { function   => 'EASIH::JMS::Samtools::get_mapped',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'get_all_unmapped' => { function   => 'EASIH::JMS::Samtools::get_unmapped',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00"},
		 
		 'identify_indel'   => { function   => 'identify_indel',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'realign_indel'    => { function   => 'realign_indel',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},

		 'realigned_merge'     => { function   => 'EASIH::JMS::Picard::merge',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=08:00:00",
					 sync       => 1},


		 'realigned_rename'     => { function   => 'rename'},
		 
		 
		 'realigned_sort'    => { function   => 'EASIH::JMS::Picard::sort',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=20000mb,walltime=08:00:00"},
		 

		 'realigned_index'   =>{ function   => 'EASIH::JMS::Samtools::index',
					hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2000mb,walltime=04:00:00"},
		 
		 'call_indels'      => { function   => 'call_indels',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},

		 'merge_indels'     => { function   => 'merge_indels',
					 sync       => 1},


		 'identify_snps'    => { function   => 'identify_snps',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=02:00:00"},
		 
		 'filter_snps'      => { function   => 'filter_snps',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500b,walltime=01:00:00"},
		 

		 'merge_vcfs'       => { function   => 'merge_vcfs',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=500mb,walltime=02:00:00", 
					 sync       => 1},
		 
		 'cluster_snps'     => { function   => 'cluster_snps',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=01:00:00"},
		 
		 'rescore_snps'     => { function   => 'rescore_snps',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=50000mb,walltime=01:00:00"},		 
		 
		 

    );
		     


our %flow = ( 'csfasta2fastq'    => 'std-aln',
	      'fastq-split'      => 'std-aln',

	      'std-aln'          => 'std-generate',
	      'std-generate'     => 'std-tag_sam',
	      'std-tag_sam'      => 'std-sam2bam',
	      'std-sam2bam'      => 'std-merge',
	      'std-merge'        => ['get_unmapped', 'get_mapped'],

	      'get_mapped'       => 'mapped_merge',
	      'get_unmapped'     => 'unmapped_bam2fq',
	      'unmapped_bam2fq'  => '2nd-aln',
	      '2nd-aln'          => '2nd-generate',
	      '2nd-generate'     => '2nd-sam2bam',
	      '2nd-sam2bam'      => 'mapped_merge',
	      'mapped_merge'     => 'mapped_sort',
	      'mapped_sort'      => 'mapped_index',

	      'mapped_index'     => ['identify_indel', 'get_all_unmapped'],
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
getopts('i:b:f:n:hlr:g:a:m:a:d:o:p:R:', \%opts);

# if ( $opts{ r} ) {
#   EASIH::JMS::restore_state($opts{r});
#   getopts('i:b:f:n:hlr:g:', \%opts);
# }


my $infile      = $opts{'i'} || $opts{'g'} || usage();
my $bam_file    = $opts{'b'} || usage();
my $reference   = $opts{'R'} || usage();
my $split       = $opts{'n'} || 10000000;
my $align_param = $opts{'a'} || ' ';
my $dbsnp       = $opts{'d'} || usage();
my $filters     = $opts{'f'} || "";
my $report      = $opts{'o'} || usage();

my $readgroup = $opts{'r'};
my $platform  = uc($opts{'p'}) || usage();


$platform = 'SOLEXA' if ( $platform eq 'ILLUMINA');
$align_param .= " -c " if ( $platform eq "SOLID");


my $bwa          = EASIH::JMS::Misc::find_program('bwa');
my $fq_split     = EASIH::JMS::Misc::find_program('fastq_split.pl');
my $samtools     = EASIH::JMS::Misc::find_program('samtools');
my $solid2fq     = EASIH::JMS::Misc::find_program('solid2fastq.pl');
my $tag_sam      = EASIH::JMS::Misc::find_program('tag_sam.pl');
my $gatk         = EASIH::JMS::Misc::find_program('gatk ');
my $sam2fq       = EASIH::JMS::Misc::find_program('sam2fastq.pl');

validate_input();

#my $solid2fq  = '/home/kb468/bin/solid2fastq.pl';
#my $bwa       = '/home/easih/bin/bwa';
#my $samtools  = '/home/easih/bin/samtools';
#my $fq_split  = '/home/kb468/bin/fastq_split.pl';


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


  # print the messages and die if critical ones.
  die join("\n", @errors) . "\n"   if ( @errors );
  print  join("\n", @info) . "\n"   if ( @info );
}





#EASIH::JMS::verbosity(10);

#EASIH::JMS::dry_run('csfasta2fastq');
#exit;
EASIH::JMS::hive('Darwin');


if ( $opts{'g'} ) {
  
  &EASIH::JMS::run('std-aln');
#  &EASIH::JMS::run('SAM2BAM');
  
}
else {
  &EASIH::JMS::run('csfasta2fastq');
}

&EASIH::JMS::store_state();

my $extra_report = "infile ==> $infile\n";
$extra_report .= "bamfile ==> $bam_file\n";
$extra_report .= "align_param ==> $align_param\n";
$extra_report .= "Binaries used..\n";
$extra_report .= `ls -l $samtools`;
$extra_report .= `ls -l $bwa` . "\n";

EASIH::JMS::mail_report('kim.brugger@easih.ac.uk', $bam_file, $extra_report);



sub csfasta2fastq {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file();
 
  my $cmd = "$solid2fq";
  $cmd .= " -n $split " if ( $split );
  $cmd .= "$infile $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}


sub fastq_split {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file();
 
  my $cmd = "$solid2fq";
  $cmd .= " -n $split " if ( $split );
  $cmd .= "$infile $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}



sub bwa_aln {
  my ($input) = @_;
  
  $input = $opts{'g'} if ($opts{'g'});
  $input =~ s/\.\z//;

  my @inputs = glob "$input";

  return if ( ! @inputs );

  foreach my $input ( @inputs ) {
    my $tmp_file = EASIH::JMS::tmp_file(".sai");
    my $cmd = "$bwa aln $align_param  -f $tmp_file $reference $input ";
    EASIH::JMS::submit_job($cmd, "$tmp_file $input");
  }
}


sub bwa_generate {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".sam");
  my $cmd = "$bwa samse -f $tmp_file $reference $input  ";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

# 
# 
# 
# Kim Brugger (04 Aug 2010)
sub bam2fq {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".fq");
  my $cmd = "$samtools view $input | $sam2fq > $tmp_file ";
  EASIH::JMS::submit_job($cmd,  $tmp_file);
}



sub bwa_aln_loose {
  my ($input) = @_;
  
  my $tmp_file = EASIH::JMS::tmp_file(".sai");
  my $cmd = "$bwa aln $align_param -e5 -t5  -f $tmp_file $reference $input ";
  EASIH::JMS::submit_job($cmd, "$tmp_file $input");
}



# 
# Adds readgroup & aligner tags to the sam-file.
# 
# Kim Brugger (22 Jul 2010)
sub sam_add_tags {
  my ($input) = @_;

  if ( ! $readgroup ) {
    $readgroup = $infile;
    $readgroup =~ s/.fastq//;
    $readgroup =~ s/.fq//;
    $readgroup =~ s/.gz//;
  }

  my $cmd = "$tag_sam -R $input -p platform -r $readgroup -a bwa -A '$align_param' ";
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
    $cmd = "$samtools view -b $bam_file $region > $tmp_file";
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


#  print "MERGE :: @inputs \n";

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
  print "snp entries: $entries\n";
  return if ( $entries == 0 );

  $filters = "--filterExpression 'DP < 20' --filterName shallow --filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' -filterName StandardFilters --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE";

#  $filters = "--filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' -filterName StandardFilters --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE";

  my $cmd = "$gatk -T VariantFiltration  -R $reference  -B variant,VCF,$input  -o $tmp_file $filters";
  EASIH::JMS::submit_job($cmd, $tmp_file);
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
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {
  
  print "Not the right usage, please look at the code\n";
  exit;

}
