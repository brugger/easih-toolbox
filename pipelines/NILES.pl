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

my $cmd = "$smalt map -f samsoft $smalt_db $first | egrep -v \\\# |  $samtools view -t $reference.fai -S - | $tag_sam -p TORRENT -r $first | $samtools view -t $reference.fai -Sb - | $samtools sort - $out; $samtools index $out.bam";
#print "$cmd\n";
system $cmd;

$cmd = "$bam_recalib -R $reference -b $out.bam -o $out.re.bam > $out.recal.csv; $samtools index $out.re.bam";
#print "$cmd\n";
system $cmd;

$cmd = "$gatk -T UnifiedGenotyper -R $reference -I $out.re.bam --max_deletion_fraction 2  | egrep -v ^INFO | ";
print "$cmd\n";
open (my $o, "> $out.vcf") || die "Could not write to file '$out.vcf': $!\n";
open (my $p, $cmd) || die "Could not open pipe: $!\n";
while(<$p>) {
  if (/\#/ ) {
    print $o $_;
    next;
  }
  print;
  my ($chr, $pos, $rest) = split("\t", $_, 3);
  next if (! $chr );

  my ($new_chr, $offset, undef) = $chr =~ /^(.*?):(\d+)-(\d+)/;
  $offset--;
  

  print $o "$new_chr\t".($pos+$offset)."\t$rest";
}
close($o);
close($p);

$cmd = "~/easih-toolbox/scripts/Variation_report.pl -s $out.vcf -O $out.var.csv -o $out.var_full.csv\n";
print "$cmd\n";
system $cmd;


# 
# 
# 
# Kim Brugger (17 Jan 2012)
sub usage {
  die "Not done yet @ARGV\n";
}




