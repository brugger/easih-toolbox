#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (01 Dec 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;


# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 1;
  if ( $DYNAMIC_LIB_PATHS ) {
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules" && ! grep /^$path\/modules/, @INC);
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    push @INC, '/home/kb468/easih-toolbox/modules/';
  }
}


#use EASIH;

my $debug = 0;
#$debug = 1;

my %barcodes;

my %opts;
getopts("f:ho:", \%opts);

usage() if ($opts{h});
my $fastq_file   = $opts{f} || usage();
my $outdir      = $opts{'o'};
$outdir = "/tmp/BC2FQ/" if ($debug);

add_barcodes();

my %bcodes;
my %fhs;

open (my $fqs, $fastq_file) || die "Could not open '$fastq_file': $!\n";
$fastq_file =~ s/\.fq//;
$fastq_file =~ s/\.gz//;

while(<$fqs>) {
  my $header     = $_;
  my $sequence   = <$fqs>;
  my $strand     = <$fqs>;
  my $quality    = <$fqs>;

  my ($bcode) = ($sequence =~ /^(.{4})/);
  next if ( ! $barcodes{$bcode});
  
  my $fout = openfile( "$fastq_file\_$barcodes{$bcode}.1");
  
  print $fout "$header$sequence$strand$quality";
  

}



# 
# 
# 
# Kim Brugger (07 Dec 2011)
sub openfile {
  my ($file_name) = @_;

  return $fhs{ $file_name } if ( $fhs{ $file_name } );

    
  my $fh;
  open ($fh, "| gzip -c > $file_name.fq.gz") || die "Could not open '$file_name': $!\n";
  $fhs{ "$file_name" }  = $fh;

  return $fhs{ $file_name }
}


# 
# 
# 
# Kim Brugger (01 Dec 2011)
sub usage {
  die "BLA BLA BLA\n";
}




# 
# 
# 
# Kim Brugger (07 Dec 2011)
sub add_barcodes {
  
  %barcodes = (  "ACGA" => "A01", 
		 "ACGC" => "A02", 
		 "AGAC" => "A03", 
		 "AGCG" => "A04", 
		 "ATCG" => "A05", 
		 "ATAC" => "A06", 
		 "CGTG" => "A07", 
		 "CTCG" => "A08", 
		 "TCTC" => "A09", 
		 "TGAC" => "A10", 
		 "TATA" => "A11", 
		 "TATC" => "A12", 
		 "TATG" => "B01", 
		 "TAGT" => "B02", 
		 "TAGA" => "B03", 
		 "TAGC" => "B04", 
		 "TACT" => "B05", 
		 "TACA" => "B06", 
		 "TACG" => "B07", 
		 "TGTA" => "B08", 
		 "TGTC" => "B09", 
		 "TGTG" => "B10", 
		 "TGAT" => "B11", 
		 "TGAG" => "B12", 
		 "TGCT" => "C01", 
		 "TGCA" => "C02", 
		 "TGCG" => "C03", 
		 "TCTA" => "C04", 
		 "TCTG" => "C05", 
		 "TCAT" => "C06", 
		 "TCAC" => "C07", 
		 "TCAG" => "C08", 
		 "TCGT" => "C09", 
		 "TCGA" => "C10", 
		 "TCGC" => "C11", 
		 "ATAT" => "C12", 
		 "ATAG" => "D01", 
		 "ATCT" => "D02", 
		 "ATCA" => "D03", 
		 "ATGT" => "D04", 
		 "ATGA" => "D05", 
		 "ATGC" => "D06", 
		 "AGTA" => "D07", 
		 "AGTC" => "D08", 
		 "AGTG" => "D09", 
		 "AGAT" => "D10", 
		 "AGAG" => "D11", 
		 "AGCT" => "D12", 
		 "AGCA" => "E01", 
		 "ACTA" => "E02", 
		 "ACTC" => "E03", 
		 "ACTG" => "E04", 
		 "ACAT" => "E05", 
		 "ACAC" => "E06", 
		 "ACAG" => "E07", 
		 "ACGT" => "E08", 
		 "GATA" => "E09", 
		 "GATC" => "E10", 
		 "GATG" => "E11", 
		 "GACT" => "E12", 
		 "GACA" => "F01", 
		 "GACG" => "F02", 
		 "GAGT" => "F03", 
		 "GAGA" => "F04", 
		 "GAGC" => "F05", 
		 "GCTA" => "F06", 
		 "GCTG" => "F07", 
		 "GCTC" => "F08", 
		 "GCAT" => "F09", 
		 "GCAG" => "F10", 
		 "GCAC" => "F11", 
		 "GCGT" => "F12", 
		 "GCGA" => "G01", 
		 "GCGC" => "G02", 
		 "CTAT" => "G03", 
		 "CTAG" => "G04", 
		 "CTAC" => "G05", 
		 "CTCT" => "G06", 
		 "CTCA" => "G07", 
		 "CTGT" => "G08", 
		 "CTGA" => "G09", 
		 "CTGC" => "G10", 
		 "CATA" => "G11", 
		 "CATC" => "G12", 
		 "CATG" => "H01", 
		 "CACT" => "H02", 
		 "CACA" => "H03", 
		 "CACG" => "H04", 
		 "CAGT" => "H05", 
		 "CAGA" => "H06", 
		 "CAGC" => "H07", 
		 "CGTA" => "H08", 
		 "CGTC" => "H09", 
		 "CGAT" => "H10", 
		 "CGAG" => "H11", 
		 "CGAC" => "H12");


}
