#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (21 Jul 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts('b:m:hg:l:L:', \%opts);

usage() if ( $opts{h});

my $bam_file  = $opts{b} || usage();
my $min_depth = $opts{m} || 200;
my $gap_file  = $opts{g};
my $low_file  = $opts{l};
my $region    = $opts{L};

my ($gout, $lout);

$gap_file = undef if ( $gap_file && $low_file &&  $gap_file eq $low_file);

open ( $gout, "> $gap_file" ) || die "Could not open file '$gap_file': $!\n" if ( $gap_file);
open ( $lout, "> $low_file" ) || die "Could not open file '$low_file': $!\n" if ( $low_file);

my $samtools  = '/usr/local/bin/samtools';
my $bam2depth = '/usr/local/bin/bam2depth';

my $seqs = names_n_lengths( $bam_file );

my (@gaps, @low_coverage);

if ( $region ) {

  my ($name, $start, $end) = $region =~ /^(\w+?):(\d+)-(\d+)/;
  die "region should be in name:start-end format not '$region'\n"
      if ( ! $name && $start && $end);

  open (my $depth, "$bam2depth $bam_file $region | ") || die "Could not open bam2depth pipeline with $bam_file $name: $!\n";
  analyse_depth( $depth, $name, $start, $end);


}
else {
  foreach my $seq ( @$seqs ) {
  
    my ($name, $length) = @$seq;
    open (my $depth, "$bam2depth $bam_file $name | ") || die "Could not open bam2depth pipeline with $bam_file $name: $!\n";
    analyse_depth( $depth, $name, 0, $length );

  }
}

print_report(\@gaps, $gout);
print_report(\@low_coverage, $lout);



# 
# 
# 
# Kim Brugger (11 Aug 2010)
sub analyse_depth {
  my ($fh, $name, $start, $end) = @_;

  my $pre_end = $start;
  my ($low_start, $low_end, $low_depth) = (-1, undef, 0);
  while (<$fh>) {
    chomp;
    my( $region, $start, $depth) = split("\t", $_);
      
    push @gaps, [$name, $pre_end + 1, $start - 1] if ( $start > $pre_end + 1 );
      
    $low_start = $start  if ( $depth < $min_depth && $low_start == -1);
    $low_end = $start    if ( $depth < $min_depth);
    $low_depth += $depth if ( $depth < $min_depth);
    
    if ( $low_start != -1 && $depth > $min_depth ) {
      push @low_coverage, [$name, $low_start, $start -1, int($low_depth/($start - $low_start))];
      $low_start = -1;
      $low_depth = 0;
    }
    
    $pre_end = $start;
  }  
  
  push @low_coverage, [$name, $low_start, $low_end, int($low_depth/($low_end - $low_start + 1))] if ( $low_start != -1 );
  
  push @gaps, [$name, $pre_end + 1, $end] if ($pre_end + 1 < $end );
}




# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub print_report {
  my ($entries, $stream ) = @_;

  $stream = *STDOUT if ( ! $stream );

  foreach my $entry ( @$entries ) {
    $$entry[3] ||= 0;
    print $stream "$$entry[0]:$$entry[1]-$$entry[2]\t$$entry[3]\n";
  }
}





# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub names_n_lengths {
  my ( $bam_file ) = @_;

  my @sequences = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);

    my ($name, $length);
    foreach my $field ( split("\t") ) {
      $name   = $1 if ( $field =~ /SN:(.*)/);
      $length = $1 if ( $field =~ /LN:(\d+)/);
    }

    push @sequences, [$name, $length] if ( $name && $length );
  }

    
  return \@sequences;
}



# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub usage {
  die "What ever, read the code...\n";
}