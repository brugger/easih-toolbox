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

  my ($bcode) = ($sequence =~ /^(.{10})/);
  next if ( !$bcode || ! $barcodes{$bcode});
  
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
  $0 =~ s/.*\///;

  print "USAGE: $0 demultiplexes a torrent run into 5 barcode plates\n";
  print "USAGE: $0 -f[astq file] -o<ut dir>\n";
  exit -1;
}




# 
# 
# 
# Kim Brugger (07 Dec 2011)
sub add_barcodes {
  
  %barcodes = ( "ACACACACAC"	=>	"p1_A1",
		"ACACACAGCA"	=>	"p1_B1",
		"ACACACATGT"	=>	"p1_C1",
		"ACACACGACG"	=>	"p1_D1",
		"ACACACGCGA"	=>	"p1_E1",
		"ACACACTAGC"	=>	"p1_F1",
		"ACACACTCTG"	=>	"p1_G1",
		"ACACACTGAT"	=>	"p1_H1",
		"ACACAGACTA"	=>	"p1_A2",
		"ACACAGAGAG"	=>	"p1_B2",
		"ACACAGCACA"	=>	"p1_C2",
		"ACACAGCGTC"	=>	"p1_D2",
		"ACACAGCTAT"	=>	"p1_E2",
		"ACACAGTCGT"	=>	"p1_F2",
		"ACACATATCG"	=>	"p1_G2",
		"ACACATCAGT"	=>	"p1_H2",
		"ACACATGATA"	=>	"p1_A3",
		"ACACATGCAG"	=>	"p1_B3",
		"ACACATGTGC"	=>	"p1_C3",
		"ACACGACACG"	=>	"p1_D3",
		"ACACGACGAC"	=>	"p1_E3",
		"ACACGACTGA"	=>	"p1_F3",
		"ACACGAGAGC"	=>	"p1_G3",
		"ACACGAGCAT"	=>	"p1_H3",
		"ACACGATATA"	=>	"p1_A4",
		"ACACGATGCT"	=>	"p1_B4",
		"ACACGCAGTC"	=>	"p1_C4",
		"ACACGCATAG"	=>	"p1_D4",
		"ACACGCGTCA"	=>	"p1_E4",
		"ACACGTACGA"	=>	"p1_F4",
		"AGACACTATA"	=>	"p1_G4",
		"AGACAGATGC"	=>	"p1_H4",
		"AGACAGCTCG"	=>	"p1_A5",
		"AGACAGTGCA"	=>	"p1_B5",
		"AGACATACGT"	=>	"p1_C5",
		"AGACATAGAC"	=>	"p1_D5",
		"AGACATCGCT"	=>	"p1_E5",
		"AGACATCTGA"	=>	"p1_F5",
		"AGACATGCTC"	=>	"p1_G5",
		"AGACGACGTG"	=>	"p1_H5",
		"AGACGACTAT"	=>	"p1_A6",
		"AGACGAGCGA"	=>	"p1_B6",
		"AGACGAGTCG"	=>	"p1_C6",
		"AGACGATAGT"	=>	"p1_D6",
		"AGACGATCAC"	=>	"p1_E6",
		"AGACGCACAT"	=>	"p1_F6",
		"AGACGCAGCG"	=>	"p1_G6",
		"AGACGCATGA"	=>	"p1_H6",
		"AGACGCGATG"	=>	"p1_A7",
		"AGACGTACTG"	=>	"p1_B7",
		"AGACGTGTAC"	=>	"p1_C7",
		"AGACTACACT"	=>	"p1_D7",
		"AGACTATATC"	=>	"p1_E7",
		"AGACTATGAG"	=>	"p1_F7",
		"AGACTCACGC"	=>	"p1_G7",
		"AGACTCAGTA"	=>	"p1_H7",
		"AGACTCGAGA"	=>	"p1_A8",
		"AGACTGAGCT"	=>	"p1_B8",
		"AGACTGATAG"	=>	"p1_C8",
		"AGACTGCATA"	=>	"p1_D8",
		"ATACAGTAGA"	=>	"p1_E8",
		"ATACAGTGTG"	=>	"p1_F8",
		"ATACATAGTA"	=>	"p1_G8",
		"ATACATCGAG"	=>	"p1_H8",
		"ATACATGTAT"	=>	"p1_A9",
		"ATACGATGTC"	=>	"p1_B9",
		"ATACGCGAGT"	=>	"p1_C9",
		"ATACGCGCAC"	=>	"p1_D9",
		"ATACGCTACA"	=>	"p1_E9",
		"ATACGCTGAG"	=>	"p1_F9",
		"ATACGTATCA"	=>	"p1_G9",
		"ATACGTCAGA"	=>	"p1_H9",
		"ATACGTGCTA"	=>	"p1_A10",
		"ATACTACATG"	=>	"p1_B10",
		"ATACTACTCA"	=>	"p1_C10",
		"ATACTATCAT"	=>	"p1_D10",
		"ATACTCATAT"	=>	"p1_E10",
		"ATACTCGATC"	=>	"p1_F10",
		"ATACTCTCTA"	=>	"p1_G10",
		"ATACTGACGA"	=>	"p1_H10",
		"ATACTGAGTC"	=>	"p1_A11",
		"ATACTGCGAT"	=>	"p1_B11",
		"ATAGACACGT"	=>	"p1_C11",
		"ATAGACTATG"	=>	"p1_D11",
		"ATAGAGATAG"	=>	"p1_E11",
		"ATAGAGCGCA"	=>	"p1_F11",
		"ATAGAGTGAC"	=>	"p1_G11",
		"ATAGATACTC"	=>	"p1_H11",
		"ATAGATCAGC"	=>	"p1_A12",
		"ATAGATGACA"	=>	"p1_B12",
		"ATCACATCTG"	=>	"p1_C12",
		"ATCACGACGA"	=>	"p1_D12",
		"ATCACGAGAC"	=>	"p1_E12",
		"ATCACGCAGC"	=>	"p1_F12",
		"ATCACGCGCA"	=>	"p1_G12",
		"ATCACGTCAT"	=>	"p1_H12",
		"ATCACTAGCG"	=>	"p2_A1",
		"ATCACTATGT"	=>	"p2_B1",
		"ATCACTCTAC"	=>	"p2_C1",
		"ATCACTGCAG"	=>	"p2_D1",
		"ATCAGACAGA"	=>	"p2_E1",
		"ATCAGACGAC"	=>	"p2_F1",
		"ATCAGAGATC"	=>	"p2_G1",
		"ATCAGAGTAG"	=>	"p2_H1",
		"ATCAGATCGC"	=>	"p2_A2",
		"ATCAGATGCA"	=>	"p2_B2",
		"ATCAGCACTG"	=>	"p2_C2",
		"ATCAGCAGAT"	=>	"p2_D2",
		"ATCAGTAGTC"	=>	"p2_E2",
		"ATCAGTCATG"	=>	"p2_F2",
		"ATCAGTCTCA"	=>	"p2_G2",
		"ATCAGTGTGC"	=>	"p2_H2",
		"ATCATACGCG"	=>	"p2_A3",
		"ATCATAGACA"	=>	"p2_B3",
		"ATCATAGCGT"	=>	"p2_C3",
		"ATCATCATAC"	=>	"p2_D3",
		"ATCATCGCTA"	=>	"p2_E3",
		"ATCATCGTCT"	=>	"p2_F3",
		"ATCATGACTC"	=>	"p2_G3",
		"ATCATGCATA"	=>	"p2_H3",
		"CACACACACG"	=>	"p2_A4",
		"CACACACGAC"	=>	"p2_B4",
		"CACACACTGA"	=>	"p2_C4",
		"CACACAGAGC"	=>	"p2_D4",
		"CACACAGCAG"	=>	"p2_E4",
		"CACACAGTCT"	=>	"p2_F4",
		"CACACATATA"	=>	"p2_G4",
		"CACACATCGT"	=>	"p2_H4",
		"CACACGACAC"	=>	"p2_A5",
		"CACACGAGCA"	=>	"p2_B5",
		"CACACGATGT"	=>	"p2_C5",
		"CACACGCATC"	=>	"p2_D5",
		"CACACGCTAG"	=>	"p2_E5",
		"CACACGTACT"	=>	"p2_F5",
		"CACACGTCTG"	=>	"p2_G5",
		"CACACTACGA"	=>	"p2_H5",
		"CACACTAGAG"	=>	"p2_A6",
		"CACACTCAGT"	=>	"p2_B6",
		"CACACTCGTA"	=>	"p2_C6",
		"CACACTGACA"	=>	"p2_D6",
		"CACACTGCTC"	=>	"p2_E6",
		"CACAGACGCA"	=>	"p2_F6",
		"CACAGACTAT"	=>	"p2_G6",
		"CACAGAGATG"	=>	"p2_H6",
		"CACAGAGCGA"	=>	"p2_A7",
		"CACAGATCAC"	=>	"p2_B7",
		"CACAGCACAG"	=>	"p2_C7",
		"CACAGCAGCT"	=>	"p2_D7",
		"CACAGCATGA"	=>	"p2_E7",
		"CACAGCGAGT"	=>	"p2_F7",
		"CGACACACAG"	=>	"p2_G7",
		"CGACACAGTC"	=>	"p2_H7",
		"CGACACATCA"	=>	"p2_A8",
		"CGACACGCGT"	=>	"p2_B8",
		"CGACACGTAC"	=>	"p2_C8",
		"CGACACTGCG"	=>	"p2_D8",
		"CGACAGAGAT"	=>	"p2_E8",
		"CGACAGCAGC"	=>	"p2_F8",
		"CGACAGCGTA"	=>	"p2_G8",
		"CGACAGTACT"	=>	"p2_H8",
		"CGACAGTCAC"	=>	"p2_A9",
		"CGACATACTA"	=>	"p2_B9",
		"CGACATCACG"	=>	"p2_C9",
		"CGACATCTAT"	=>	"p2_D9",
		"CGACATGAGA"	=>	"p2_E9",
		"CGACGACATA"	=>	"p2_F9",
		"CGACGACGCT"	=>	"p2_G9",
		"CGACGACTGC"	=>	"p2_H9",
		"CGACGAGCTC"	=>	"p2_A10",
		"CGACGATACG"	=>	"p2_B10",
		"CGACGCGACT"	=>	"p2_C10",
		"CGACGCTAGA"	=>	"p2_D10",
		"CGACGCTCTG"	=>	"p2_E10",
		"CGACGCTGAC"	=>	"p2_F10",
		"CGACGTACAC"	=>	"p2_G10",
		"CGACGTAGCA"	=>	"p2_H10",
		"CGACGTATGT"	=>	"p2_A11",
		"CGACGTCGAG"	=>	"p2_B11",
		"CGACTACGTC"	=>	"p2_C11",
		"CGACTAGAGC"	=>	"p2_D11",
		"CTACACACGA"	=>	"p2_E11",
		"CTACACGAGC"	=>	"p2_F11",
		"CTACACGCTG"	=>	"p2_G11",
		"CTACACTCAT"	=>	"p2_H11",
		"CTACACTGTA"	=>	"p2_A12",
		"CTACAGAGCG"	=>	"p2_B12",
		"CTACAGATAC"	=>	"p2_C12",
		"CTACAGCATG"	=>	"p2_D12",
		"CTACAGCTCA"	=>	"p2_E12",
		"CTACATCTGC"	=>	"p2_F12",
		"CTACATGACT"	=>	"p2_G12",
		"CTACATGCAC"	=>	"p2_H12",
		"CTACGACAGT"	=>	"p3_A1",
		"CTACGACTAG"	=>	"p3_B1",
		"CTACGAGATG"	=>	"p3_C1",
		"CTACGAGTGA"	=>	"p3_D1",
		"CTACGATCTA"	=>	"p3_E1",
		"CTACGATGAT"	=>	"p3_F1",
		"CTACGCACTC"	=>	"p3_G1",
		"CTACGCAGCT"	=>	"p3_H1",
		"CTACGCGTAT"	=>	"p3_A2",
		"CTACGTAGTG"	=>	"p3_B2",
		"CTACGTCTCT"	=>	"p3_C2",
		"CTACGTGCGT"	=>	"p3_D2",
		"CTACTAGTAC"	=>	"p3_E2",
		"CTACTATAGA"	=>	"p3_F2",
		"CTACTCAGAG"	=>	"p3_G2",
		"CTACTCATGC"	=>	"p3_H2",
		"CTACTCGACA"	=>	"p3_A3",
		"CTACTCTATG"	=>	"p3_B3",
		"CTCACATGCG"	=>	"p3_C3",
		"CTCACGCTCT"	=>	"p3_D3",
		"CTCACTACTG"	=>	"p3_E3",
		"CTCAGAGTGT"	=>	"p3_F3",
		"CTCAGTACGT"	=>	"p3_G3",
		"CTCAGTAGCA"	=>	"p3_H3",
		"CTCAGTCACT"	=>	"p3_A4",
		"CTCAGTGAGA"	=>	"p3_B4",
		"CTCATACTAC"	=>	"p3_C4",
		"CTCATAGCTG"	=>	"p3_D4",
		"CTCATATGTA"	=>	"p3_E4",
		"CTCATCGCAT"	=>	"p3_F4",
		"CTCATCTATC"	=>	"p3_G4",
		"CTCATCTCGA"	=>	"p3_H4",
		"CTCATGATCA"	=>	"p3_A5",
		"CTCATGTCAC"	=>	"p3_B5",
		"CTCGAGATCT"	=>	"p3_C5",
		"CTCGAGTGTG"	=>	"p3_D5",
		"CTCGATCTGT"	=>	"p3_E5",
		"CTCGCAGTGA"	=>	"p3_F5",
		"CTCGCGCACA"	=>	"p3_G5",
		"CTCGCGTGCT"	=>	"p3_H5",
		"CTCGCTACGC"	=>	"p3_A6",
		"CTCGTATAGT"	=>	"p3_B6",
		"CTCGTCAGCG"	=>	"p3_C6",
		"CTCGTCGATA"	=>	"p3_D6",
		"CTCGTGACGA"	=>	"p3_E6",
		"CTCGTGCATG"	=>	"p3_F6",
		"CTCTACAGAT"	=>	"p3_G6",
		"CTCTACTCGT"	=>	"p3_H6",
		"TACACTCGCG"	=>	"p3_A7",
		"TACACTGTAC"	=>	"p3_B7",
		"TACAGAGTGC"	=>	"p3_C7",
		"TACAGATATC"	=>	"p3_D7",
		"TACAGATGAG"	=>	"p3_E7",
		"TACAGCAGTG"	=>	"p3_F7",
		"TACAGCTCGA"	=>	"p3_G7",
		"TACAGTACAC"	=>	"p3_H7",
		"TACAGTCGTC"	=>	"p3_A8",
		"TACAGTCTCT"	=>	"p3_B8",
		"TACAGTGACG"	=>	"p3_C8",
		"TACAGTGCTA"	=>	"p3_D8",
		"TACATAGAGT"	=>	"p3_E8",
		"TACATAGTAG"	=>	"p3_F8",
		"TACATATACA"	=>	"p3_G8",
		"TACATATCGC"	=>	"p3_H8",
		"TACATCATCA"	=>	"p3_A9",
		"TACATGACGA"	=>	"p3_B9",
		"TACATGAGAC"	=>	"p3_C9",
		"TACATGCTAT"	=>	"p3_D9",
		"TACATGTATG"	=>	"p3_E9",
		"TACGAGCAGC"	=>	"p3_F9",
		"TACGAGTCTC"	=>	"p3_G9",
		"TACGATACTG"	=>	"p3_H9",
		"TACGATATCT"	=>	"p3_A10",
		"TACGATCGAG"	=>	"p3_B10",
		"TACGCACATA"	=>	"p3_C10",
		"TACGCACTCT"	=>	"p3_D10",
		"TACGCAGCGT"	=>	"p3_E10",
		"TACGCATAGC"	=>	"p3_F10",
		"TCACATCTAG"	=>	"p3_G10",
		"TCACATGTCA"	=>	"p3_H10",
		"TCACGATCTC"	=>	"p3_A11",
		"TCACGCGTGT"	=>	"p3_B11",
		"TCACGCTGCA"	=>	"p3_C11",
		"TCACGTAGCG"	=>	"p3_D11",
		"TCACTAGTCG"	=>	"p3_E11",
		"TCACTATAGC"	=>	"p3_F11",
		"TCACTATGTG"	=>	"p3_G11",
		"TCACTCACAG"	=>	"p3_H11",
		"TCACTCATCT"	=>	"p3_A12",
		"TCACTCGATG"	=>	"p3_B12",
		"TCACTCGCGC"	=>	"p3_C12",
		"TCACTGATAC"	=>	"p3_D12",
		"TCACTGTACA"	=>	"p3_E12",
		"TCAGACTATC"	=>	"p3_F12",
		"TCAGAGAGTC"	=>	"p3_G12",
		"TCAGAGATCG"	=>	"p3_H12",
		"TCAGAGCGAT"	=>	"p4_A1",
		"TCAGATAGAG"	=>	"p4_B1",
		"TCAGATATGT"	=>	"p4_C1",
		"TCAGATGCTA"	=>	"p4_D1",
		"TCAGCATCTG"	=>	"p4_E1",
		"TCAGCGCAGC"	=>	"p4_F1",
		"TCAGCGTGCG"	=>	"p4_G1",
		"TCAGCTCATG"	=>	"p4_H1",
		"TCAGCTCTCA"	=>	"p4_A2",
		"TCAGCTGAGT"	=>	"p4_B2",
		"TCAGCTGTAC"	=>	"p4_C2",
		"TCAGTACTAT"	=>	"p4_D2",
		"TGACACATAT"	=>	"p4_E2",
		"TGACACGACA"	=>	"p4_F2",
		"TGACACTCTC"	=>	"p4_G2",
		"TGACAGACTG"	=>	"p4_H2",
		"TGACAGCGAG"	=>	"p4_A3",
		"TGACAGTCGA"	=>	"p4_B3",
		"TGACATGATG"	=>	"p4_C3",
		"TGACATGTGT"	=>	"p4_D3",
		"TGACGCACTA"	=>	"p4_E3",
		"TGACGCGCAG"	=>	"p4_F3",
		"TGACGTATAG"	=>	"p4_G3",
		"TGACGTCACA"	=>	"p4_H3",
		"TGACGTGAGC"	=>	"p4_A4",
		"TGACTACAGA"	=>	"p4_B4",
		"TGACTACGAT"	=>	"p4_C4",
		"TGACTATCTA"	=>	"p4_D4",
		"TGACTCTACG"	=>	"p4_E4",
		"TGACTCTCAT"	=>	"p4_F4",
		"TGACTGCTCA"	=>	"p4_G4",
		"TGACTGTAGT"	=>	"p4_H4",
		"TGAGACGCAC"	=>	"p4_A5",
		"TGAGACTAGA"	=>	"p4_B5",
		"TGAGATCGTC"	=>	"p4_C5",
		"TGAGCACAGT"	=>	"p4_D5",
		"TGAGCAGTAG"	=>	"p4_E5",
		"TGAGCGACGC"	=>	"p4_F5",
		"TGAGCGTCAG"	=>	"p4_G5",
		"TGAGCTAGAC"	=>	"p4_H5",
		"TGAGCTCTGC"	=>	"p4_A6",
		"TGAGTACGTA"	=>	"p4_B6",
		"TGCACAGCTA"	=>	"p4_C6",
		"TGCACGACGT"	=>	"p4_D6",
		"TGCACGCAGA"	=>	"p4_E6",
		"TGCACGTCAC"	=>	"p4_F6",
		"TGCACTATGA"	=>	"p4_G6",
		"TGCACTCACT"	=>	"p4_H6",
		"TGCACTCTAG"	=>	"p4_A7",
		"TGCACTGCAT"	=>	"p4_B7",
		"TGCAGCGTAT"	=>	"p4_C7",
		"TGCAGCTGCT"	=>	"p4_D7",
		"TGCAGTGTCA"	=>	"p4_E7",
		"TGCATACGCA"	=>	"p4_F7",
		"TGCATAGACG"	=>	"p4_G7",
		"TGCATATGTG"	=>	"p4_H7",
		"TGCATCGCGT"	=>	"p4_A8",
		"TGCATCTCTC"	=>	"p4_B8",
		"TGCATGAGCT"	=>	"p4_C8",
		"TGCGACGATC"	=>	"p4_D8",
		"TGCGACGTGT"	=>	"p4_E8",
		"TGCGACTGCG"	=>	"p4_F8",
		"TGCGAGACAG"	=>	"p4_G8",
		"TGCGAGCGTA"	=>	"p4_H8",
		"TGCGATGAGA"	=>	"p4_A9",
		"TGCGCACGTG"	=>	"p4_B9",
		"TGCGCAGACA"	=>	"p4_C9",
		"TGCGCATCTC"	=>	"p4_D9",
		"TGCGCGATCT"	=>	"p4_E9",
		"TGCGCTAGCA"	=>	"p4_F9",
		"TGCGCTCATC"	=>	"p4_G9",
		"TGCGTAGAGC"	=>	"p4_H9",
		"ATGACATGTC"	=>	"p4_A10",
		"ATGACGAGTG"	=>	"p4_B10",
		"ATGACGATAT"	=>	"p4_C10",
		"ATGACGTAGA"	=>	"p4_D10",
		"ATGACTATCA"	=>	"p4_E10",
		"ATGACTCGAT"	=>	"p4_F10",
		"ATGAGACTGC"	=>	"p4_G10",
		"ATGAGATACT"	=>	"p4_H10",
		"ATGAGATGAG"	=>	"p4_A11",
		"ATGAGCACGT"	=>	"p4_B11",
		"ATGAGCGATA"	=>	"p4_C11",
		"ATGAGCTAGC"	=>	"p4_D11",
		"ATGAGTACAC"	=>	"p4_E11",
		"ATGAGTCTAG"	=>	"p4_F11",
		"ATGAGTGACG"	=>	"p4_G11",
		"ATGATACAGT"	=>	"p4_H11",
		"ATGATAGCAC"	=>	"p4_A12",
		"ATGATAGTCG"	=>	"p4_B12",
		"ATGATATCTA"	=>	"p4_C12",
		"ATGATCAGCG"	=>	"p4_D12",
		"ATGATCTCAG"	=>	"p4_E12",
		"ATGATGATGA"	=>	"p4_F12",
		"ATGATGCGAC"	=>	"p4_G12",
		"ATGATGTATC"	=>	"p4_H12",
		"ATGCACACAG"	=>	"p5_A1",
		"ATGCACATCA"	=>	"p5_B1",
		"ATGCACGAGA"	=>	"p5_C1",
		"ATGCACTACT"	=>	"p5_D1",
		"ATGCACTCGC"	=>	"p5_E1",
		"ATGCAGACGT"	=>	"p5_F1",
		"ATGCAGAGAC"	=>	"p5_G1",
		"ATGCAGCAGC"	=>	"p5_H1",
		"ATGCAGCGCG"	=>	"p5_A2",
		"ATGCAGTCTA"	=>	"p5_B2",
		"ATGCATCGTC"	=>	"p5_C2",
		"ATGCATGATG"	=>	"p5_D2",
		"ATGCGACGAT"	=>	"p5_E2",
		"ATGCGACTCG"	=>	"p5_F2",
		"ATGCGAGACA"	=>	"p5_G2",
		"ATGCGAGCGT"	=>	"p5_H2",
		"CTGACACATA"	=>	"p5_A3",
		"CTGACACGAG"	=>	"p5_B3",
		"CTGACACTGT"	=>	"p5_C3",
		"CTGACAGACG"	=>	"p5_D3",
		"CTGACAGCGC"	=>	"p5_E3",
		"CTGACGCGTC"	=>	"p5_F3",
		"CTGACGTATG"	=>	"p5_G3",
		"CTGACTACAT"	=>	"p5_H3",
		"CTGACTCAGC"	=>	"p5_A4",
		"CTGACTCGCA"	=>	"p5_B4",
		"CTGACTGCTA"	=>	"p5_C4",
		"CTGACTGTAG"	=>	"p5_D4",
		"CTGAGAGCAT"	=>	"p5_E4",
		"CTGAGATATC"	=>	"p5_F4",
		"CTGAGCATAC"	=>	"p5_G4",
		"CTGAGCGCTG"	=>	"p5_H4",
		"CTGAGCGTGA"	=>	"p5_A5",
		"CTGAGCTGCA"	=>	"p5_B5",
		"CTGAGTATCT"	=>	"p5_C5",
		"CTGAGTCGTG"	=>	"p5_D5",
		"CTGATATCGT"	=>	"p5_E5",
		"CTGATATGAC"	=>	"p5_F5",
		"CTGATCATGT"	=>	"p5_G5",
		"CTGATCGACT"	=>	"p5_H5",
		"CTGATCTGTG"	=>	"p5_A6",
		"CTGATGACGC"	=>	"p5_B6",
		"CTGATGAGCT"	=>	"p5_C6",
		"CTGCACGTAG"	=>	"p5_D6",
		"CTGCACTGAC"	=>	"p5_E6",
		"CTGCAGCTAT"	=>	"p5_F6",
		"CTGCAGTGCA"	=>	"p5_G6",
		"CTGCATAGAT"	=>	"p5_H6",
		"CTGCATCATA"	=>	"p5_A7",
		"CTGCATCTCG"	=>	"p5_B7",
		"CTGCATGTGA"	=>	"p5_C7",
		"CTGCGACGCA"	=>	"p5_D7",
		"CTGCGCATCG"	=>	"p5_E7",
		"CTGCGCGATC"	=>	"p5_F7",
		"CTGCGCTAGT"	=>	"p5_G7",
		"CTGCGCTCAG"	=>	"p5_H7",
		"TGTACACGAT"	=>	"p5_A8",
		"TGTACAGTGT"	=>	"p5_B8",
		"TGTACATCTG"	=>	"p5_C8",
		"TGTACATGCA"	=>	"p5_D8",
		"TGTACGAGTA"	=>	"p5_E8",
		"TGTACTACTC"	=>	"p5_F8",
		"TGTACTCATG"	=>	"p5_G8",
		"TGTAGAGTAG"	=>	"p5_H8",
		"TGTAGCACAG"	=>	"p5_A9",
		"TGTAGCAGTC"	=>	"p5_B9",
		"TGTAGCGAGT"	=>	"p5_C9",
		"TGTAGTACGA"	=>	"p5_D9",
		"TGTAGTCTGC"	=>	"p5_E9",
		"TGTATAGATC"	=>	"p5_F9",
		"TGTATCATGA"	=>	"p5_G9",
		"TGTATCTGAG"	=>	"p5_H9",
		"TGTATGCATA"	=>	"p5_A10",
		"TGTATGCTCT"	=>	"p5_B10",
		"TGTCACAGTA"	=>	"p5_C10",
		"TGTCACGCTG"	=>	"p5_D10",
		"TGTCACGTGA"	=>	"p5_E10",
		"TGTCAGATCG"	=>	"p5_F10",
		"TGTCATGACT"	=>	"p5_G10",
		"TGTCATGTAG"	=>	"p5_H10",
		"TGTCGACTCT"	=>	"p5_A11",
		"TGTCGCGTCG"	=>	"p5_B11",
		"TGTCGCTACT"	=>	"p5_C11",
		"TGTCGTAGTG"	=>	"p5_D11",
		"TGTCTACATG"	=>	"p5_E11",
		"TGTCTAGCGT"	=>	"p5_F11",
		"TGTCTCGATA"	=>	"p5_G11",
		"TGTCTGACGA"	=>	"p5_H11",
		"TGTCTGCGAC"	=>	"p5_A12",
		"TGTCTGTCAG"	=>	"p5_B12",
		"TGTGACACGT"	=>	"p5_C12",
		"TGTGAGCTGA"	=>	"p5_D12",
		"TGTGATAGCT"	=>	"p5_E12",
		"TGTGATGTCA"	=>	"p5_F12",
		"TGTGCATACT"	=>	"p5_G12",
		"TGTGCGCTAC"	=>	"p5_H12",
      );
}
