package EASIH::Barcodes::illumina;
# 
# 
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use vars qw ( $bc_length
              %barcodes);



$bc_length = 6;


# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub bc_length {
  return $bc_length;
}

%barcodes = (
  "ATCACG"  =>   1,
  "CGATGT"  =>   2, 
  "TTAGGC"  =>   3,
  "TGACCA"  =>   4,
  "ACAGTG"  =>   5,
  "GCCAAT"  =>   6,
  "CAGATC"  =>   7,
  "ACTTGA"  =>   8,
  "GATCAG"  =>   9,
  "TAGCTT"  =>  10,
  "GGCTAC"  =>  11,
  "CTTGTA"  =>  12,
  "AGTCAA"  =>  13,
  "AGTTCC"  =>  14,
  "ATGTCA"  =>  15,
  "CCGTCC"  =>  16,
  "GTAGAG"  =>  17,
  "GTCCGC"  =>  18,
  "GTGAAA"  =>  19,
  "GTGGCC"  =>  20,
  "GTTTCG"  =>  21,
  "CGTACG"  =>  22,
  "GAGTGG"  =>  23,
  "GGTAGC"  =>  24,
  "ACTGAT"  =>  25,
  "ATGAGC"  =>  26,
  "ATTCCT"  =>  27,
  "CAAAAG"  =>  28,
  "CAACTA"  =>  29,
  "CACCGG"  =>  30,
  "CACGAT"  =>  31,
  "CACTCA"  =>  32,
  "CAGGCG"  =>  33,
  "CATGGC"  =>  34,
  "CATTTT"  =>  35,
  "CCAACA"  =>  36,
  "CGGAAT"  =>  37,
  "CTAGCT"  =>  38,
  "CTATAC"  =>  39,
  "CTCAGA"  =>  40,
  "GACGAC"  =>  41,
  "TAATCG"  =>  42,
  "TACAGC"  =>  43,
  "TATAAT"  =>  44,
  "TCATTC"  =>  45,
  "TCCCGA"  =>  46,
  "TCGAAG"  =>  47,
  "TCGGCA"  =>  48,
  # NEB barcodes::
  "CGTGAT"  => "NEBNext1",
  "ACATCG"  => "NEBNext2",
  "GCCTAA"  => "NEBNext3",
  "TGGTCA"  => "NEBNext4",
  "CACTGT"  => "NEBNext5",
  "ATTGGC"  => "NEBNext6",
  "GATCTG"  => "NEBNext7",
  "TCAAGT"  => "NEBNext8",
  "CTGATC"  => "NEBNext9",
  "AAGCTA"  => "NEBNext10",
  "GTAGCC"  => "NEBNext11",
  "TACAAG"  => "NEBNext12");

# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub barcodes {
  return %barcodes;
  
}



1;



