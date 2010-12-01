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



our %analysis = ('identify_snps'    => { function   => 'identify_snps',
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
		     


our %flow = ( 'identify_snps'    => 'filter_snps',
	      'filter_snps'      => 'merge_vcfs',
	      'merge_vcfs'       => 'cluster_snps',
	      'cluster_snps'     => 'rescore_snps'
	      );

#EASIH::JMS::no_store();
#EASIH::JMS::print_flow('fastq-split');

my %opts;
getopts('b:R:d:f:o:hH:L:p:M:S:', \%opts);

usage() if ( $opts{h});
my $hard_reset    = $opts{'H'};
my $soft_reset    = $opts{'S'};

if ( $soft_reset ) {
  print "Doing a soft reset/restart\n";
  &EASIH::JMS::reset($soft_reset);
  getopts('b:R:d:f:o:hH:L:p:M:S:', \%opts);
}
elsif ( $hard_reset ) {
  &EASIH::JMS::hard_reset($hard_reset);
  getopts('b:R:d:f:o:hH:L:p:M:S:', \%opts);
}


my $bam_file    = $opts{'b'} || usage();
my $reference   = $opts{'R'} || usage();
my $dbsnp       = $opts{'d'} || usage();
my $filter      = $opts{'f'} || "exon";
my $log         = $opts{'L'};
open (*STDOUT, ">> $log") || die "Could not open '$log': $!\n" if ( $log );
my $min_depth   = $opts{'M'}     || 0;
my $report      = $opts{'o'} || usage();

my $platform    = uc($opts{'p'});
$platform = 'SOLEXA'      if ( $platform eq 'ILLUMINA');


my $samtools     = EASIH::JMS::Misc::find_program('samtools');
my $gatk         = EASIH::JMS::Misc::find_program('gatk');

validate_input();

#EASIH::JMS::verbosity(10);
EASIH::JMS::backend('Darwin');
#EASIH::JMS::hive('Kluster');
EASIH::JMS::max_retry(0);

if ($filter eq "wgs" ) {
  $min_depth ||= 20;

  $filter  = "--filterExpression 'DP < $min_depth' --filterName shallow ";
#  $filter .= "-clusterWindowSize 10 ";
  $filter .= "--filterExpression 'AB > 0.75 && DP > 40 || DP > 100 || MQ0 > 40 || SB > -0.10'   --filterName StandardFilters ";
  $filter .= "--filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)'  --filterName HARD_TO_VALIDATE ";
  
}
elsif ( $filter eq "wgs-low" ) {
  $min_depth ||= 5;

  $filter  = "--filterExpression 'DP < $min_depth' --filterName shallow ";
  $filter .= "--filterExpression 'MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1'  --filterName HARD_TO_VALIDATE";
}
elsif ( $filter eq "exon" ) {
  $min_depth ||= 20;

  $filter = "--filterExpression 'DP < $min_depth' --filterName shallow ";
  $filter .= " --filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' -filterName StandardFilters ";
  $filter .= " --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE";  
}
elsif ( $filter eq "exon-low" ) {
  $min_depth ||= 5;

  $filter = "--filterExpression 'DP < $min_depth' --filterName shallow ";
  $filter .= " --filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' -filterName StandardFilters ";
  $filter .= " --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE";  
}


&EASIH::JMS::run('identify_snps');

&EASIH::JMS::store_state();

my $extra_report .= "bamfile ==> $bam_file\n";
$extra_report .= "snp_file ==> $report.snps\n";

$extra_report .= "Binaries used..\n";
$extra_report .= `ls -l $samtools`;

EASIH::JMS::mail_report('kim.brugger@easih.ac.uk', $bam_file, $extra_report);

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
  return if ( $entries == 0 );

  my $cmd = "$gatk -T VariantFiltration  -R $reference  -B variant,VCF,$input  -o $tmp_file $filter";
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
  
  push @errors, "'$dbsnp' does not exists\n" if (! -e $dbsnp);
  push @errors, "'$dbsnp' does end with .rod as expected\n" if ($dbsnp !~ /.rod\z/);

  push @errors, "Platform must be either SOLEXA or SOLID not '$platform'" if ( $platform && $platform ne "SOLEXA" && $platform ne 'SOLID');

  push @errors, "$filter should be one of the following: wgs,wgs-low,exon,exon-low\n"
      if ($filter ne "wgs" && $filter ne "wgs-low" && $filter ne "exon" && $filter ne "exon-low");

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
  print "USAGE: $0 -b[am file] -R [eference genome] -d[bsnp rod] -o[ut prefix] -p[latform: illumina or solid\n\n";

  print "extra flags: -f[ilter: wgs,wgs-low,exon,exon-low. Default= exon] \n";
  print "extra flags: -H[ard reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "extra flags: -L[og file, default is STDOUT]\n";  
  print "extra flags: -M[in depth for snps, defaults: normal=20 low=5]\n";
  print "extra flags: -S[oft reset/restart of a crashed/failed run, needs a freeze file]\n";
  print "\n";
  print "easih-pipeline: " . &EASIH::JMS::version() . "\n";
  exit;

}
