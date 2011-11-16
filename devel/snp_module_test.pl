#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (08 Dec 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::SNPs;

EASIH::SNPs->New('dbsnp_134');
print "Hellow world, let look at some snps...\n";

my $snp = EASIH::SNPs::fetch_snp(1, 11370591);

print Dumper( $snp );

my $cm = EASIH::SNPs::CM(1, 11370691);
$cm = EASIH::SNPs::CM(1, 11377114);

print "At 1:11370691 $cm\n";
