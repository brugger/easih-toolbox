#!/usr/bin/perl 
# 
# Generic script for snp analysis.
# 
# 
# Kim Brugger (07 Jan 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my @homozygous;
my @heterozygous;
my $max_snp_freq    = 10;
my $min_sample_size = 30;
my $background;
my $KGPilot         = 0;
my $novel           = 0;

GetOptions ("Homozygous=s"    => \@homozygous,
	    "heterozygous=s"  => \@heterozygous,
	    "max_snp_freq=i"  => \$max_snp_freq,
	    "min_sample_size" => \$min_sample_size,
	    "background=s"    => \$background,
	    "KGPilot"         => \$KGPilot,
	    "novel"           => \$novel,
    );

@homozygous   = split(/,/,join(',',@homozygous));
@heterozygous = split(/,/,join(',',@heterozygous));


my (%SNPs,%all_SNPs);

my (%field_names,$region, @field_names);

my ($inputs, $file1) = (0);


foreach my $file (@homozygous, @heterozygous) {

  open( my $in, $file) || die "Could not open '$file': $!\n";

  $inputs++;
  my $escape = 10;

  while(<$in>) {
    next if (/#/);
    chomp;

    my @f = split("\t");
    
    if ( ! %field_names &&  $f[0] eq "position") {
      @field_names = @f;
      $file1 = $file;
      for(my $i=0; $i< @f; $i++ ) {
	
	$field_names{ $f[$i] } = $i;

	if ($f[$i] eq 'region' ) { 
 	  $region = $i;
	  $region = 13;
	}
      }
      next;
    }
    else { 
      my $effect = $f[$region];

      next if ($SNPs{$file}{$f[0]});

      
      if ( $effect eq "ESSENTIAL_SPLICE_SITE" ||
	   $effect eq "NON_SYNONYMOUS_CODING" ||
	   $effect eq "SPLICE_SITE" ||
	   $effect eq "STOP_GAINED" ||
	   $effect eq "STOP_LOST" ) {

	$all_SNPs{$f[0]}++;
	
	@{$SNPs{$file}{$f[0]}} = @f;

      }
    }
  }

}

my $goodies = 0;
foreach my $key (keys %all_SNPs) {
  if ( $all_SNPs{ $key} != $inputs ) {
#    print "$all_SNPs{ $key} != $inputs\n";
    delete $all_SNPs{ $key};
  }
  else {
    $goodies++;
  }
}
print "#$goodies harmful SNPs are shared between all the $inputs sample(s)\n";

#
# Remove snps due to frequency.
#
$goodies = 0;
foreach my $key (keys %all_SNPs) {

  if ( ${$SNPs{$file1}{$key}}[ $field_names{ 'dbsnp freq' }] ) {
    my ( $samples, $freq) = split('/', ${$SNPs{$file1}{$key}}[ $field_names{ 'dbsnp freq' }]);
    
    $freq *= 100;

    if ($min_sample_size < $samples &&
	$freq            < $max_snp_freq ) {
      $goodies++;
    }
    else {
      delete $all_SNPs{ $key};
    }
  }
  else {
    $goodies++;
  }

}

print "#$goodies harmful SNPs left after removing SNPs with a frequency < $max_snp_freq\n";


if ( $novel ) {
  $goodies = 0;
  foreach my $key (keys %all_SNPs) {

    if ( ${$SNPs{$file1}{$key}}[ $field_names{ 'external ref' }] &&
	 ${$SNPs{$file1}{$key}}[ $field_names{ 'external ref' }] ne "") {
      delete $all_SNPs{ $key};
    }
    else {
      $goodies++;
    }
  }

  print "#$goodies harmful SNPs left after removing known SNPs from snpdb\n";
}

if ( $KGPilot ) { 
  $goodies = 0;
  foreach my $key (keys %all_SNPs) {

    if ( ${$SNPs{$file1}{$key}}[ $field_names{ 'external ref' }] &&
	 ${$SNPs{$file1}{$key}}[ $field_names{ 'dbsnp flags' }]  &&
	 (${$SNPs{$file1}{$key}}[ $field_names{ 'dbsnp flags' }] eq "" ||
	  ${$SNPs{$file1}{$key}}[ $field_names{ 'dbsnp flags' }] !~ /KGPilot/ )) {
      $goodies++;
    }
    else {
      delete $all_SNPs{ $key};
    }
  }
  print "#$goodies harmful SNPs left after removing all 1kg SNPs\n";
}

if ( $background ) {

  $goodies = 0;
  open (my $in, $background) || die "Could not open '$background': $!\n";
  while(<$in>) {
    
    chomp;
    
    if ( $all_SNPs{ $_ }) {
      delete $all_SNPs{ $_ };
    }
  }
  
  print "#" .(keys %all_SNPs) ." harmful SNPs left after removing from EASIH-SNP-DB\n";
}


$goodies = 0;
foreach my $file ( @homozygous ) {

  foreach my $key (keys %all_SNPs) {
   if ( $SNPs{$file}{$key} &&
	${$SNPs{$file}{$key}}[$field_names{ 'type' }] !~ /homozygous/) {
     delete $all_SNPs{ $key };
   }
   else {
     $goodies++;
   }
  }
}


$goodies = 0;

foreach my $file ( @heterozygous ) {

  foreach my $key (keys %all_SNPs) {
   if ( $SNPs{$file}{$key} &&
	${$SNPs{$file}{$key}}[$field_names{ 'type' }] !~ /heterozygous/) {
     delete $all_SNPs{ $key };
   }
   else {
     $goodies++;
   }
  }
}

print "#" .(keys %all_SNPs) ." harmful SNPs left after using the homozygous/heterozygous logic\n";

#exit;

print join("\t", @field_names) . "\n";

foreach my $key (sort keys %all_SNPs) {
  print join("\t", @{$SNPs{$file1}{$key}}) . "\n";
}

