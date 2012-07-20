# 
# 
# 
# 
# Kim Brugger (07 Feb 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use Test::Simple tests => 19;

use lib '/home/kb468/easih-toolbox/modules/';

use EASIH::Sample;


my $good_sample = EASIH::Sample::validate_sample("AA00001");
ok(! $good_sample, 'Failing invalid name: AA00001' );

$good_sample = EASIH::Sample::validate_sample("A00001");
ok(! $good_sample, 'Failing invalid name: A00001' );

$good_sample = EASIH::Sample::validate_sample("0000001");
ok(! $good_sample, 'Failing invalid name: 000001' );

$good_sample = EASIH::Sample::validate_sample("a000001");
ok(! $good_sample, 'Failing invalid name: a000001' );

$good_sample = EASIH::Sample::validate_sample("A000001_");
ok(! $good_sample, 'Failing invalid name: A000001_' );

$good_sample = EASIH::Sample::validate_sample("A000001");
ok($good_sample, 'Validating: A000001' );

$good_sample = EASIH::Sample::validate_sample("A000001a");
ok($good_sample, 'Validating: A000001a' );

$good_sample = EASIH::Sample::validate_sample("A0000001");
ok($good_sample, 'Validating: A0000001' );

$good_sample = EASIH::Sample::validate_sample("A0000001a");
ok($good_sample, 'Validating: A0000001a' );

my $project = EASIH::Sample::extract_project("A00001");
ok(! defined $project, 'Cannot extract project from invalid sample name' );

$project = EASIH::Sample::extract_project("A0100001a");
ok( $project eq "A01", 'Extracted project name' );

my ($sample, $version, $postfix) = EASIH::Sample::sample_n_version("A0100001a.1.fq.gz");
ok( $sample eq "A0100001a" && $version == 0 && $postfix eq ".1.fq.gz", 'Correctly split A0100001a.1.fq.gz');

($sample, $version, $postfix) = EASIH::Sample::sample_n_version("A0100001a_1.1.fq.gz");
ok( $sample eq "A0100001a" && $version == 1 && $postfix eq ".1.fq.gz", 'Correctly split A0100001a_1.1.fq.gz');

($sample, $version, $postfix) = EASIH::Sample::sample_n_version("A0100001a_1");
ok( $sample eq "A0100001a" && $version == 1 && $postfix eq "", 'Correctly split A0100001a_1');

my $next_sample_name = EASIH::Sample::next_sample_name('A010001');
ok( $next_sample_name eq "A010001", 'Correctly named next sample for A010001, no prior filenames');

$next_sample_name = EASIH::Sample::next_sample_name('A010001', ['A010001']);
ok( $next_sample_name eq "A010001_1", 'Correctly named next sample for A010001, one prior filenames');

$next_sample_name = EASIH::Sample::next_sample_name('A010001', ['A010001_4', 'A010001']);
ok( $next_sample_name eq "A010001_5", 'Correctly named next sample for A010001, two prior filenames');

$next_sample_name = EASIH::Sample::next_sample_name('A010001', ['A010002']);
ok( $next_sample_name eq "A010001", 'Correctly named next sample for A010001, with one wrong prior filenames');

my ($sample_file, $error) = EASIH::Sample::sample2outfilename( $next_sample_name );
ok( $sample_file eq "/data/A01/raw/A010001", 'Correctly generated file path from sample name');
