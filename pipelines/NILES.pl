#!/usr/bin/perl 
# 
# pipeline for few genes analysis, mainly for clinical use.
# 
# 
# Kim Brugger (17 Jan 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

BEGIN {
  use vars qw/$path/; 
  $path = $0;
  if ($path =~ /.*\//) {
    $path =~ s/(.*\/).*/$1/;
  }
  else {
    $path = "./";
  }
  push @INC, $path;
}

use Getopt::Std;

#use lib '/home/cjp64/git/easih-pipeline/modules';
use lib '/home/kb468/easih-pipeline/modules';
use lib '/home/kb468/easih-toolbox/modules';

use EASIH::Misc;

my $opts = '1:R:s:p:';
my %opts;
getopts($opts, \%opts);

#usage() if ( $opts{h});


my $first         = $opts{'1'}     || usage();
#$first            = [split(",", $first)];
die "Can only handle single ends reads for now, sorry\n" if ( $opts{'2'} );
my $platform      = uc($opts{'p'}) || usage();
$platform = 'ILLUMINA'      if ( $platform eq 'ILLUMINA');
my $reference    = $opts{'R'}     || usage();
my $smalt_db     = $opts{'s'}     || usage();

my $samtools    = EASIH::Misc::find_program('samtools');
my $smalt       = EASIH::Misc::find_program('smalt');
my $tag_sam     = EASIH::Misc::find_program('tag_sam.pl');
my $bam_recalib = EASIH::Misc::find_program('bam_recalibrate.pl');
my $gatk        = EASIH::Misc::find_program('gatk_1.3-14');


my $smalt_ref = $reference;
$smalt_ref =~ s/^(.*)\..*/$1/;

my $out = $first;

#$out =~ s/^([A-Z]\d{6,7}).*/$1/;
$out =~ s/.1.fq//;
$out =~ s/.gz//;

#my $cmd = "$smalt map -f samsoft $smalt_db $first | egrep -v \\\# |  $samtools view -t $reference.fai -S - | $tag_sam -p TORRENT -r $first | $samtools view -t $reference.fai -Sb - | $samtools sort - $out; $samtools index $out.bam";
my $cmd = "$smalt map -f samsoft $smalt_db $first | egrep -v \\\# | $tag_sam -p TORRENT -r $first | ";
#print "$cmd\n";
#system $cmd;

my %regions;

open( my $bam_out, " | $samtools view -t $reference.fai -Sb - | $samtools sort - $out" ) || die "Could not open stream: $!\n";
open( my $bam_in,  $cmd ) || die "Could not open stream: $!\n";
while(<$bam_in>) {
  if (/^\@/) {
    print $bam_out $_;
  }
  else {
    my @F = split("\t");
    if ( $F[2] =~ /(.*?):(\d+)-(\d+)/) {
      my ($chr, $start, $end) = ($1,$2,$3);
      $regions{ $F[2] }++;
      $F[3] += $start - 1;
      $F[2] = $chr;
      print $bam_out join("\t", @F);
    }
  }
}
close($bam_out);
close($bam_in);
system "$samtools index $out.bam";

$cmd = "$bam_recalib -R $reference -b $out.bam -o $out.re.bam > $out.recal.csv; $samtools index $out.re.bam";
print "$cmd\n";
system $cmd;


open (my $var_out, " > $out.intervals") || die "Could not open file '$out.intervals': $!\n";
foreach my $region (keys %regions) {
  print $var_out "$region\n";
}
close( $var_out );

$cmd = "$gatk -T UnifiedGenotyper -R $reference -I $out.re.bam --max_deletion_fraction 2 -o $out.re.vcf -L $out.intervals ";
print "$cmd\n";
system $cmd;

$cmd = "~/easih-toolbox/scripts/Variation_report.pl -s $out.re.vcf -O $out.var.csv -o $out.var_full.csv\n";
print "$cmd\n";
system $cmd;
 

# 
# 
# 
# Kim Brugger (17 Jan 2012)
sub usage {
  die "Not done yet @ARGV\n";
}




