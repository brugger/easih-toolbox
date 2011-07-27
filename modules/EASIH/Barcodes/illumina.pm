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
  "ATCACG" => 1,
  "CGATGT" => 2, 
  "TTAGGC" =>  3,
  "TGACCA" =>  4,
  "ACAGTG" =>  5,
  "GCCAAT" =>  6,
  "CAGATC" =>  7,
  "ACTTGA" =>  8,
  "GATCAG" =>  9,
  "TAGCTT" =>  10,
  "GGCTAC" =>  11,
  "CTTGTA" =>  12);



# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub barcodes {
  return %barcodes;
  
}



1;



