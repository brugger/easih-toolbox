#!/usr/bin/perl 
# 
# rmdup script that will work better for fragments than what picard is currently doing
# 
# 
# Kim Brugger (03 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts('b:h', \%opts);
my $bam_file = $opts{'b'} || Usage();
#my $sam_file = $opts{'s'} || Usage();

my $samtools    = find_program('samtools');


my (%seen, $pre_pos);
open (my $bam, "$samtools view -h  $bam_file | ") || die "Could not open 'stream': $!\n";
while(<$bam>) {

  if (/^\@/) {
    print;
    next;
  }

  my @F = split("\t");
  my ($id, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, @opts) = @F;


  if (1 || $flags & 0x0004) {

    
    if ( !$pre_pos || $pos != $pre_pos) {
      %seen = ();
      
      $seen{ $sequence }{ $cigar }++;
      $pre_pos = $pos;
    }
    else {

      if ( $seen{ $sequence }{ $cigar } ) {

	$F[1] += 0x400;
      }
      else {
	$seen{ $sequence }{ $cigar }++;
      }
    }

  }
    print join("\t",@F);

  
}


# 
# 
# 
# Kim Brugger (13 Jul 2010)
sub find_program {
  my ($program) = @_;
  
  my $username = scalar getpwuid $<;
  
  my @paths = ("/home/$username/bin/",
	       "./",
	       "/usr/local/bin");
  
  foreach my $path ( @paths ) {
    return "$path/$program" if ( -e "$path/$program" );
  }

  my $location = `which $program`;
  chomp( $location);
  
  return $location if ( $location );
  
  die "Could not find '$program'\n";
}


# 
# 
# 
# Kim Brugger (05 Nov 2010)
sub Usage {
  $0 =~ s/.*\///;
  die "USAGE: $0 -b<am file>\n";

  # tests :::: odd SNP reporting:  10:74879852-74879852
  # large indel: 10:111800742
}
