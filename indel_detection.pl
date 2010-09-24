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



our %analysis = ('call_indels'      => { function   => 'call_indels',
					 hpc_param  => "-NEP-fqs -l nodes=1:ppn=1,mem=2500b,walltime=12:00:00"},

		 'merge_indels'     => { function   => 'merge_indels',
					 sync       => 1},
		 );
		     


our %flow = ( 'call_indels'      => 'merge_indels',
	      );

#EASIH::JMS::no_store();
#EASIH::JMS::print_flow('fastq-split');

my %opts;
getopts('b:R:o:p:h', \%opts);

# if ( $opts{ r} ) {
#   EASIH::JMS::restore_state($opts{r});
#   getopts('i:b:f:n:hlr:g:', \%opts);
# }


my $bam_file    = $opts{'b'} || usage();
my $reference   = $opts{'R'} || usage();
my $filters     = $opts{'f'} || "default";
my $report      = $opts{'o'} || usage();

my $readgroup   = $opts{'r'} || $report;

my $samtools     = EASIH::JMS::Misc::find_program('samtools');
my $gatk         = EASIH::JMS::Misc::find_program('gatk ');

validate_input();

#EASIH::JMS::verbosity(10);
EASIH::JMS::backend('Darwin');
#EASIH::JMS::backend('Kluster');
EASIH::JMS::max_retry(0);


&EASIH::JMS::run('call_indels');

&EASIH::JMS::store_state();

my $extra_report .= "bamfile ==> $bam_file\n";
$extra_report .= "indel_file ==> $report.indel\n";

$extra_report .= "Binaries used..\n";
$extra_report .= `ls -l $samtools`;

EASIH::JMS::mail_report('kim.brugger@easih.ac.uk', $bam_file, $extra_report);

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
  print "USAGE: $0 -b [am file] -R [eference genome]  -o[ut prefix] \n";
  exit;

}
