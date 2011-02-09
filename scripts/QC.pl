#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (28 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

use strict;
use Getopt::Long;

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::QC;

my %opts;
getopts('b:c:f:q:hs:', \%opts);

my $fastq_file    = $opts{f};
my $csfastq_file  = $opts{c};
my $qual_file     = $opts{q};
my $bam_file      = $opts{b};
my $sample_size   = $opts{s} || 20;

#sample limit MB
EASIH::QC::sample_size($sample_size);

if ( $fastq_file ) {
  my $QC = EASIH::QC::fastQC( $fastq_file, 1 );
  $fastq_file =~ s/(.*?.[12])\..*/$1/ ||  $fastq_file =~ s/(.*?)\..*/$1/; 
  EASIH::QC::make_plots($QC, $fastq_file);
  print "Predicted mappability : $$QC{mappability} %\n" if ( $$QC{mappability} );
}

if ( $csfastq_file ) {
  my $QC = EASIH::QC::csfastaQC( $csfastq_file);
  $csfastq_file =~ s/(.*?)\..*/$1/;
  EASIH::QC::make_plots($QC, $csfastq_file);
}

if ( $qual_file ) {
  my $QC = EASIH::QC::qualQC( $qual_file );
  $qual_file =~ s/(.*?)\..*/$1/;
  EASIH::QC::make_plots($QC, $qual_file);
  print "Predicted mappability : $$QC{mappability} %\n" if ( $$QC{mappability} );
  print "+Q30: $$QC{Q30}\n" if ($$QC{Q30});
}

if ( $bam_file ) {
  my ($QC1, $QC2) = EASIH::QC::bamQC( $bam_file );
  $bam_file =~ s/(.*?)\..*/$1/;
  
  EASIH::QC::make_plots($QC1, "$bam_file.1");
  print "Predicted mappability : $$QC1{mappability} %\n" if ( $$QC1{mappability} );

  EASIH::QC::make_plots($QC2, "$bam_file.2");
  print "Predicted mappability : $$QC2{mappability} %\n" if ( $$QC2{mappability} );
}



__END__

# 
# 
# 
# Kim Brugger (18 Aug 2010)
sub usage {
  
  $0 =~ s/.*\///;
  print "There are three running modes, one where the seq is sampled and aligned and the other from a pre aligned file, finally just qv and base dist.\n";
  print "USAGE RAW: $0 -1 [first fq file] -2 [second fq] -R[eference bwa formatted] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  print "USAGE RAW: $0 -b[am file] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  print "USAGE RAW: $0 -1 [first fq file] -2 [second fq] -n[o mapping] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  exit;

}
