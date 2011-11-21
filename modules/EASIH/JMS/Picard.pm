package EASIH::JMS::Picard;

use strict;
use warnings;
use Data::Dumper;

use EASIH::JMS;
use EASIH::JMS::Misc;

my $picard   = EASIH::JMS::Misc::find_program('picard');

sub merge { 
  my (@inputs) = @_;

  print Dumper( \@inputs);

  @inputs = @{$inputs[0]} if ( @inputs == 1 && ref($inputs[0]) eq "ARRAY" );

  print Dumper( \@inputs);

  my $tmp_file = EASIH::JMS::tmp_file(".merged.bam");

#  print "MERGE :: @inputs \n";

  if (@inputs == 1 ) {
    EASIH::JMS::submit_system_job("mv @inputs $tmp_file", $tmp_file);
  }
  else {

    # remove empty files as they crash the merging step.
    my @non_empty_files;
    foreach my $input ( @inputs ) {
      push @non_empty_files, $input if ( ! -z $input );
    }

    my $username = scalar getpwuid $<;

    my $cmd = "$picard -T MergeSamFiles USE_THREADING=true O= $tmp_file  I= " . join(" I= ", @non_empty_files) . " VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/";
    EASIH::JMS::submit_job($cmd, $tmp_file);
  }
}



sub sort { 
  my ($input) = @_;

  my $username = scalar getpwuid $<;

  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd = "$picard -T SortSam  I= $input O= $tmp_file SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub fixmate_n_sort  { 
  my ($input) = @_;

  my $username = scalar getpwuid $<;
  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd = "$picard -T FixMateInformation  I= $input O= $tmp_file SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}

sub fixmate  { 
  my ($input) = @_;

  my $username = scalar getpwuid $<;
  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $cmd = "$picard -T FixMateInformation  I= $input O= $tmp_file  VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/";
  EASIH::JMS::submit_job($cmd, $tmp_file);
}


# 
# 
# 
# Kim Brugger (20 Sep 2010)
sub mark_duplicates {

  my ($input) = @_;

  my $username = scalar getpwuid $<;
  my $tmp_file = EASIH::JMS::tmp_file(".bam");
  my $metrix_file = EASIH::JMS::tmp_file(".mtx");
  my $cmd = "$picard -T MarkDuplicates  I= $input O= $tmp_file  M= $metrix_file VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/$username/scratch/tmp/ MAX_RECORDS_IN_RAM=500000";
  EASIH::JMS::submit_job($cmd, $tmp_file);

}



# 
# 
# 
# Kim Brugger (08 Nov 2010)
sub version {

  my $version = `$picard -v`;
  chomp($version);
  return $version;
}




1;
