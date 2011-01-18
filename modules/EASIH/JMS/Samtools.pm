package EASIH::JMS::Samtools;
#
# 
# 
# 
# 
# Kim Brugger (13 Jul 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;


use EASIH::JMS;
use EASIH::JMS::Misc;

my $samtools    = EASIH::JMS::Misc::find_program('samtools');
my $fix_sorted  = EASIH::JMS::Misc::find_program('bam_fix_sorted_flag.pl');


sub sam2bam {
  my ($input) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd = "$samtools view -b -S $input > $tmp_file ";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub merge { 
  my (@inputs) = @_;

  my $tmp_file = EASIH::JMS::tmp_file(".merged.bam");

  print "MERGE :: @inputs \n";

  if (@inputs == 1 ) {
    EASIH::JMS::submit_system_job("mv @inputs $tmp_file", $tmp_file);
  }
  else {
    my $cmd = "$samtools merge $tmp_file @inputs ";
    print "$cmd \n";
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}

sub sort {
  my ($input) = @_;

  my $tmp_file  = EASIH::JMS::tmp_file(".sorted");
  my $cmd = "$samtools sort -m 2048000000 $input $tmp_file; $fix_sorted $tmp_file";
    
  EASIH::JMS::submit_job($cmd, "$tmp_file.bam");
}


sub index {
  my ($input) = @_;

  my $cmd = "$samtools index $input";
  EASIH::JMS::submit_job($cmd, $input);
}		     


sub get_mapped {
  my ($input) = @_;
  
  my $tmp_file = EASIH::JMS::tmp_file(".mapped.bam");
  my $cmd = "$samtools view -bF4 $input > $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub get_unmapped {
  my ($input) = @_;
  
  my $tmp_file = EASIH::JMS::tmp_file(".unmapped.bam");
  my $cmd = "$samtools view -bf4 $input > $tmp_file";

  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub version {
  open( my $pipe, "$samtools 2>&1 | head -n 3 |") || die "Could not open samtools pipe: $!\n";
  <$pipe>;
  <$pipe>;
  my $version = <$pipe>;
  chomp( $version);

  return( $version);
}

1;
