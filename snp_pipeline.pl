#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use Getopt::Std;

use lib '/home/kb468/projects/easih-flow/modules';
use lib '/home/kb468/easih-pipeline/modules';
use EASIH::JMS;
use EASIH::JMS::Misc;


our %analysis = ('identify_snps' => { function   => 'identify_snps',
				      hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=00:12:00"},
		 
		 'filter_snps'   => { function   => 'filter_snps',
				      hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=00:12:00"},
		 

		 'merge_vcfs'    => { function   => 'merge_vcfs',
				      hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=00:20:00", 
				      sync       => 1},
		 
		 'cluster_snps'  => { function   => 'cluster_snps',
				      hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500mb,walltime=00:20:00"},
		 
		 'rescore_snps'  => { function   => 'rescore_snps',
				      hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=50000mb,walltime=00:30:00"},
    );
		     

our %flow = ( 'identify_snps' => 'filter_snps',
	      'filter_snps'   => 'merge_vcfs',
	      'merge_vcfs'    => 'cluster_snps',
	      'cluster_snps'  => 'rescore_snps');

my %opts;
getopts('b:R:p:f:d:o:hr:', \%opts);

if ( $opts{ r} ) {
  EASIH::JMS::restore_state($opts{r});
  getopts('b:R:p:f:d:o:hr:', \%opts);
}


my $bam_file   = $opts{'b'} || usage();
my $dbsnp      = $opts{'d'} || usage();
my $filters    = $opts{'f'} || usage();
my $platform   = $opts{'p'};
my $report     = $opts{'o'} || usage();
my $reference  = $opts{'R'} || usage();

my $resources  = " -resources ~/src/gatk/Sting/R/  -Rscript /usr/bin/Rscript ";

my $samtools  = EASIH::JMS::Misc::find_program('samtools');
my $gatk      = EASIH::JMS::Misc::find_program('gatk ');

EASIH::JMS::hive('Darwin');
#EASIH::JMS::hive('Kluster');
EASIH::JMS::max_retry(0);

&EASIH::JMS::run('identify_snps');

&EASIH::JMS::store_state();


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

  my $tmp_file = EASIH::JMS::tmp_file(".filtered.vcf");
  $filters = "--filterExpression 'DP < 20' --filterName shallow --filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' -filterName StandardFilters --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE";

  $filters = "--filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' -filterName StandardFilters --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE";

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
  my $cmd = "$gatk -T GenerateVariantClusters  -R $reference -B input,VCF,$input --DBSNP $dbsnp -an QD -an SB -an HaplotypeScore -an HRun -clusterFile $tmp_file  ";
  EASIH::JMS::submit_job($cmd, "$tmp_file -B input,VCF,$input");
}		
		
		

sub rescore_snps {
  my ($input) = @_;

  $report =~ s/.vcf\z//;
  my $cmd = "$gatk -T VariantRecalibrator -R $reference --DBSNP $dbsnp -clusterFile $input -output  $report --target_titv $resources ";
  EASIH::JMS::submit_job($cmd, $report);
}		







# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {
  
  print "Not the right usage, please look at the code\n";
  exit;

}

