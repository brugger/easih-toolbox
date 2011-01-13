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

print "Hellow world, let look at some snps...\n";

print Dumper(EASIH::SNPs::phylop_score_GRCh37(10, 60007));
print Dumper(EASIH::SNPs::phylop_score_GRCh37(10, 6000));
print Dumper(EASIH::SNPs::phylop_score_GRCh37(10, 60002));
print Dumper(EASIH::SNPs::phylop_score_GRCh37(10, 60009));
print Dumper(EASIH::SNPs::phylop_score_GRCh37(10, 60037));

exit;

print Dumper(EASIH::SNPs::fetch_snp_hg18(1, 536205));
print Dumper(EASIH::SNPs::fetch_snp_GRCh37(1, 536205));

print Dumper(EASIH::SNPs::fetch_flags('rs75786306'));

exit;


print Dumper(EASIH::SNPs::population_stats('rs1048488', 'ASW'));
print Dumper(EASIH::SNPs::population_stats('rs1048488'));



print Dumper(EASIH::SNPs::fetch_rs('rs6690870'));


