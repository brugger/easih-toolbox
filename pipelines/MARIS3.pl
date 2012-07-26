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

use lib '/scratch/kb468//BRC/easih-toolbox/modules';
use lib '/scratch/kb468//BRC/easih-pipeline/modules';

use EASIH::Pipeline;
use EASIH::Pipeline::Misc;


my $opts = '1:2:d:De:f:hH:I:lL:M:n:No:p:Q:Pr:R:sS:vV';
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

  $opts{Q} =~ s/\..*//;

  $opts{'1'} = join(",", sort(glob("$opts{Q}*.1.fq"), glob("$opts{Q}*.1.fq.gz")));
  $opts{'2'} = join(",", sort(glob("$opts{Q}*.2.fq"), glob("$opts{Q}*.2.fq.gz")));

  $opts{'L'} = "$opts{Q}.log";
  $opts{'o'} = "$opts{Q}";
  $opts{'l'} = 1;
  
  $opts{'r'} ||= '/scratch/easih/GATK_bundle/b37/';

}  

my $first          = $opts{'1'}     || usage();
$first             = [split(",", $first)];
my $second         = $opts{'2'}     || usage();
$second            = [split(",", $second)];
my $email          = $opts{'e'}     || "$username\@cam.ac.uk";
my $loose_mapping  = $opts{'l'}     || 0;
my $log            = $opts{'L'} || "$opts{o}.log" || undef;
our $report        = $opts{'o'}     || usage();
my $platform       = 'ILLUMINA';
my $reference_dir = $opts{'r'}     || usage();
my $small_capture = $opts{'s'} || 0;

my $reference   = glob("$reference_dir/*.fasta");
my ($dbsnp)     = grep(!/excluding/, glob("$reference_dir/dbsnp_*.vcf"));
my ($hapmap)    = glob("$reference_dir/hapmap_*.sites.vcf");
my ($omni)      = glob("$reference_dir/*omni*.sites.vcf");

my $bam_file      = "$report.bam";
my $vcf_file      = "$report.vcf";
my $host_cpus     = nr_of_cpus();

#$host_cpus = 1;

my $freeze_file = "$opts{'o'}.maris";
system "mv $freeze_file $freeze_file.backup"  if ( -e $freeze_file );

EASIH::Pipeline::freeze_file($freeze_file);



#EASIH::Pipeline::verbosity(100) if ( $opts{v});

open (*STDOUT, ">> $log") || die "Could not open '$log': $!\n" if ( $log );

# set platform specific bwa aln parameters
my $align_param .= " -q 15 ";
# and loose mapping
$align_param    .= " -e5 "     if ( $loose_mapping);

my $bwa             = EASIH::Pipeline::Misc::find_program('bwa_0.6.1-tpx');
#$bwa = "/home/easih/bin/bwa_0.6.2";
my $samtools        = EASIH::Pipeline::Misc::find_program('samtools');
my $gatk            = EASIH::Pipeline::Misc::find_program('gatk');
my $picard          = EASIH::Pipeline::Misc::find_program('picard');

$gatk = "/home/kb468/bin/gatk";


#validate_input();

#print "$bwa\n$samtools\n$gatk\n$picard\n";


#EASIH::Pipeline::verbosity(10);
#EASIH::Pipeline::backend('Darwin');
#EASIH::Pipeline::backend('MPIexec');
EASIH::Pipeline::backend('Local');
EASIH::Pipeline::max_jobs( $host_cpus );
EASIH::Pipeline::max_retry(0);

EASIH::Pipeline::add_start_step('bwa_aln');
EASIH::Pipeline::add_step('bwa_aln', 'bwa_sampe', );
EASIH::Pipeline::add_step('bwa_sampe', 'bam_sort');
EASIH::Pipeline::add_merge_step('bam_sort', 'bam_merge');
EASIH::Pipeline::add_step('bam_merge', 'mark_dups');
EASIH::Pipeline::add_step('mark_dups', 'realign_targets');
EASIH::Pipeline::add_step('realign_targets', 'bam_realign');
EASIH::Pipeline::add_merge_step('bam_realign', 'bam_merge2' );
EASIH::Pipeline::add_step('bam_merge2', 'bam_index2');
EASIH::Pipeline::add_step('bam_index2', 'count_covariates');
EASIH::Pipeline::add_step('count_covariates', 'table_recalibration');
EASIH::Pipeline::add_step('table_recalibration', 'bam_index3');
EASIH::Pipeline::add_step('bam_index3', 'UnifiedGenotyper_par');
EASIH::Pipeline::add_merge_step('UnifiedGenotyper_par', 'CombineVariants1', 'CombineVariants');
EASIH::Pipeline::add_step('CombineVariants1', 'VariantAnnotator');
EASIH::Pipeline::add_step('UnifiedGenotyper', 'VariantAnnotator');
EASIH::Pipeline::add_step('VariantAnnotator', 'SelectSNPVariants');
EASIH::Pipeline::add_step('SelectSNPVariants', 'VariantRecalibrator');
EASIH::Pipeline::add_step('VariantRecalibrator', 'ApplyRecalibration');
EASIH::Pipeline::add_merge_step('ApplyRecalibration', 'CombineVariants');
EASIH::Pipeline::add_step('VariantAnnotator', 'SelectIndelVariants');
EASIH::Pipeline::add_step('SelectIndelVariants', 'VariantFiltration');
EASIH::Pipeline::add_merge_step('VariantFiltration', 'CombineVariants');
EASIH::Pipeline::add_step('CombineVariants', 'run_stats');
EASIH::Pipeline::add_step('run_stats', 'finished');

#EASIH::Pipeline::print_flow();

my $run_id = get_run_id();

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
    
    my $tmp_sai1  = EASIH::Pipeline::tmp_file(".sai");
    my $tmp_sai2  = EASIH::Pipeline::tmp_file(".sai");
    
    my $cmd = "$bwa aln -t $host_cpus $align_param -f $tmp_sai1 $reference $first_file ; ";
    $cmd   .= "$bwa aln -t $host_cpus $align_param -f $tmp_sai2 $reference $second_file;";
    
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
  $cmd = "$bwa sampe -t $host_cpus -P $reference $$input{first_sai} $$input{second_sai} $$input{first_fq} $$input{second_fq} 2> tmp.err | $samtools view -Sb - -o $tmp_file";


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

  my $cmd = "$picard -T AddOrReplaceReadGroups.jar I=$input O=$tmp_file SORT_ORDER=coordinate CN=EASIH ID=$readgroup PL=$platform LB=$readgroup PU=$run_id  SM=$sample VALIDATION_STRINGENCY=SILENT  TMP_DIR=/home/$username/scratch/tmp/";

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}



# 
# 
# 
# Kim Brugger (22 Jul 2012)
sub get_run_id {
  my ($input) = @_;
  $input ||= $$first[0];

  my $line1 = ` zcat -f $input  | head -n 1 `;
  chomp $line1;
  my ($id, $rest) = split(/:/, $line1,2);
  $id =~ s/^\@//;
  
  return $id;
  
}



# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_merge {
  my (@inputs) = @_;


  @inputs = @{$inputs[0]} if ( @inputs == 1 && ref($inputs[0]) eq "ARRAY" );


  my $tmp_file = EASIH::Pipeline::tmp_file(".merged.bam");

  if (@inputs == 1 ) {
#    print "cd; mv @inputs $tmp_file, $tmp_file\n";
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
# Kim Brugger (25 Jun 2012)
sub bam_merge2 {
  my (@inputs) = @_;


  @inputs = @{$inputs[0]} if ( @inputs == 1 && ref($inputs[0]) eq "ARRAY" );


  my $tmp_file = EASIH::Pipeline::tmp_file(".merged.bam");

  if (@inputs == 1 ) {
    
#    print "cd; mv @inputs $tmp_file, $tmp_file\n";
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
sub mark_dups {
  my ($input) = @_;

  my $username = scalar getpwuid $<;
  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");
  my $metrix_file = EASIH::Pipeline::tmp_file(".mtx");
  my $cmd = "$picard -T MarkDuplicates  I= $input O= $tmp_file  M= $metrix_file VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/ CREATE_INDEX=true";
  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}



# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub realign_targets {
  my ($input) = @_;

  my @names = ();
  open(my $spipe, "$samtools view -H $input | ") || die "Could not open '$input': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }

  print "$input @names\n";

  foreach my $name ( @names ) {
    my $tmp_file = EASIH::Pipeline::tmp_file(".intervals");
    my $cmd = "$gatk -T RealignerTargetCreator -R $reference -o $tmp_file  -L $name -I $input";
    EASIH::Pipeline::submit_job($cmd, "$tmp_file $name $input");
  }
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub bam_realign {
  my ($input) = @_;

  my ($interval_file, $region, $tmp_bam_file) = split(" ", $input);

  my $tmp_file = EASIH::Pipeline::tmp_file(".bam");
  my $cmd;
  # If the interval file is empty the realigner ignores the region and produces an empty bamfile...
  if (  -z $interval_file ) {
    $cmd = "$samtools view -b $tmp_bam_file $region > $tmp_file";
  }
  else {
    $cmd = "$gatk -T IndelRealigner -targetIntervals $interval_file -L $region -o $tmp_file -R $reference -I $tmp_bam_file ";
  }

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}

# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub count_covariates {
  my ($input) = @_;
  
  my $tmp_file = EASIH::Pipeline::tmp_file("_recal.csv.");

  my $cmd = "$gatk -T CountCovariates -R $reference -I $input -cov ReadGroupCovariate  -cov QualityScoreCovariate  -cov CycleCovariate  -cov DinucCovariate  -recalFile $tmp_file ";
  $cmd .= " -knownSites $dbsnp " if ( $dbsnp);

  EASIH::Pipeline::submit_job($cmd, "$input $tmp_file");
}

# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub table_recalibration {
  my ($input) = @_;
  my ($tmp_bam, $recal) = split(" ", $input);

  my $cmd = "$gatk -T TableRecalibration  -R $reference -I $tmp_bam -recalFile $recal -baq RECALCULATE -o $bam_file"; 

  EASIH::Pipeline::submit_job($cmd, $bam_file);
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub UnifiedGenotyper {
  my ($input) = @_;

  my $tmp_file = EASIH::Pipeline::tmp_file(".vcf");
  my $cmd = "$gatk -T UnifiedGenotyper -nt $host_cpus -R $reference -I $input -glm BOTH -G Standard -A AlleleBalance -stand_call_conf 30.0 -stand_emit_conf 10.0 -dcov 1000 -baq CALCULATE_AS_NECESSARY -o $tmp_file";
  $cmd .= " --dbsnp $dbsnp ";
  

  EASIH::Pipeline::submit_job($cmd, $tmp_file);
}


# 
# 
# 
# Kim Brugger (06 Jul 2012), 
sub UnifiedGenotyper_par {
  my ($input) = @_;

  my @names = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      push @names, $1 if ( $field =~ /SN:(.*)/);
    }
  }


  foreach my $name ( @names ) {
    my $tmp_file = EASIH::Pipeline::tmp_file(".vcf");

    my $cmd = "$gatk -T UnifiedGenotyper -R $reference -I $input -glm BOTH -G Standard -A AlleleBalance -stand_call_conf 30.0 -stand_emit_conf 10.0 -dcov 1000 -baq CALCULATE_AS_NECESSARY -o $tmp_file -L $name ";
    $cmd .= " --dbsnp $dbsnp ";


    EASIH::Pipeline::submit_job($cmd, $tmp_file);

  }
}




# 
# 
# 
# Kim Brugger (26 Jun 2012)
sub VariantAnnotator {
  my ($input) = @_;


  my $tmp_vcf = EASIH::Pipeline::tmp_file("_annot.vcf");
  my $cmd = "$gatk -T VariantAnnotator -R $reference -I $bam_file -o $tmp_vcf -V $input -L $input --dbsnp $dbsnp";

  EASIH::Pipeline::submit_job($cmd, $tmp_vcf);
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub SelectSNPVariants {
  my ($input) = @_;

  my $tmp_snp_file = EASIH::Pipeline::tmp_file(".vcf");

  EASIH::Pipeline::submit_job("$gatk -T SelectVariants -R $reference -V $input -selectType SNP -o $tmp_snp_file", $tmp_snp_file);
}


# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub SelectIndelVariants {
  my ($input) = @_;

  my $tmp_indel_file = EASIH::Pipeline::tmp_file(".vcf");

  EASIH::Pipeline::submit_job("$gatk -T SelectVariants -R $reference -V $input -selectType INDEL -o $tmp_indel_file", $tmp_indel_file);
  
}




# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub VariantRecalibrator {
  my ($input) = @_;

  my $tmp_recal    = EASIH::Pipeline::tmp_file(".recal");
  my $tmp_tranches = EASIH::Pipeline::tmp_file(".tranches");
  my $tmp_r_script = EASIH::Pipeline::tmp_file(".plots.R");

  my $cmd = "$gatk -T VariantRecalibrator -R $reference -input $input ";
  $cmd .= "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap " if ( $hapmap );
  $cmd .= "-resource:omni,known=false,training=true,truth=false,prior=12.0 $omni " if ( $omni );
  $cmd .= "-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbsnp " if ( $dbsnp);
  $cmd .= "-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ ";
  $cmd .= " -recalFile $tmp_recal ";
  $cmd .= " -tranchesFile $tmp_tranches ";
  $cmd .= " -rscriptFile $tmp_r_script ";

  if ( $small_capture ) {
    $cmd .= "--minNumBadVariants 100 ";
    $cmd .= "  -mG 2 ";
  }

  EASIH::Pipeline::submit_job($cmd, "$input $tmp_recal $tmp_tranches $tmp_r_script");

}




# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub ApplyRecalibration {
  my ($input ) = @_;

  my ($vcf, $tmp_recal, $tmp_tranches, $tmp_r_script) = split(" ", $input);

  my $tmp_recal_vcf = EASIH::Pipeline::tmp_file("_recal.vcf");

  my $cmd = "$gatk -T ApplyRecalibration -R $reference -input $vcf --ts_filter_level 99.0 -tranchesFile $tmp_tranches -recalFile $tmp_recal -o $tmp_recal_vcf";

  EASIH::Pipeline::submit_job($cmd, $tmp_recal_vcf);
}



sub CombineVariants {
  my ($inputs) = @_;

  if ( @$inputs == 1 ) {
    EASIH::Pipeline::submit_system_job("mv @$inputs $vcf_file", $vcf_file);
  }
  else {
    my $cmd = "$gatk -T CombineVariants -R $reference -o $vcf_file ";
    my $count = 1;
    foreach my $input ( @$inputs ) {
      $cmd .= " -V:variant $input ";
      $count++;
    }
    EASIH::Pipeline::submit_job($cmd, $vcf_file);
  }
}





# 
# 
# 
# Kim Brugger (25 Jun 2012)
sub VariantFiltration {
  my ($input) = @_;

  my $tmp_vcf = EASIH::Pipeline::tmp_file("_indel.vcf");
  my $cmd = "$gatk -T VariantFiltration  -R $reference -V $input -filter 'QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8  || FS > 200.0' -filterName GATKStandard -o $tmp_vcf ";

  EASIH::Pipeline::submit_job($cmd, $tmp_vcf);

}



sub run_stats {
  my ($input) = @_;
  
  EASIH::Pipeline::submit_job("$samtools flagstat $report.bam > $report.flagstat");
  EASIH::Pipeline::submit_job("md5sum $bam_file > $bam_file.md5 ");
  EASIH::Pipeline::submit_job("md5sum $report.bam.bai > $report.bam.bai.md5 ");
  EASIH::Pipeline::submit_job("md5sum $report.vcf > $report.vcf.md5 ");
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
