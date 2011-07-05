#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (23 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment. Needs to be prior to the use of EASIH* modules.
BEGIN {
  my $path = $0;
  if ($path =~ /.*\//) {
    $path =~ s/(.*)\/.*/$1/;
    push @INC, "$path/modules" if ( -e "$path/modules");
    $path =~ s/(.*)\/.*/$1/;
    push @INC, "$path/modules" if ( -e "$path/modules");
    
  }
  else {
    push @INC, "../modules" if ( -e "../modules");
    push @INC, "./modules" if ( -e "./modules");
  }
}


use EASIH;
use EASIH::QC;
use EASIH::QC::db;


my $sid = EASIH::QC::db::add_sample('Z980001', 'Z98');
my $pid = EASIH::QC::db::fetch_project_id('Z98');
my $fid = EASIH::QC::db::add_file('Z980001.1.fq', 'Z980001', 'Z98');
print "$fid --> $pid --> $sid\n";

# add data functions
#EASIH::QC::db::add_qvs_histogram($fid, [[1,10],[20,200], [35,250]]);
#EASIH::QC::db::add_qvs_boxplot($fid, [[1,0,1,2,3,4,5],[2,0,1,2,3,4,5], [3,0,1,2,3,4,5]]);
#EASIH::QC::db::add_adaptors($fid, [[3, 34.4], [4, 99.4]]);
#EASIH::QC::db::add_duplicate_seqs($fid, [["AA", 34.4, "ill"], ["BBB", 0.994, "Solid"]]);
#EASIH::QC::db::add_duplicates($fid, [[3, 34.4], [4, 99.4]]);
#EASIH::QC::db::add_gc_distribution($fid, [[5, 34.4], [6, 99.4]]);
#EASIH::QC::db::add_basedists($fid, [[1, 25,25,25,25,0], [2, 27,23,20,30,0]]);


# and fetch the data again
my @qvs = EASIH::QC::db::fetch_qvs_boxplot($fid);
print Dumper( \@qvs);

my @qvs_hist = EASIH::QC::db::fetch_qvs_histogram($fid);
print Dumper( \@qvs_hist);

my @adapt = EASIH::QC::db::fetch_adaptors($fid);
print Dumper( \@adapt );

my @dup_seqs = EASIH::QC::db::fetch_duplicate_seqs($fid);
print Dumper(\@dup_seqs);

my @dups = EASIH::QC::db::fetch_duplicates($fid);
print Dumper(\@dups );

my @gc_dist = EASIH::QC::db::fetch_gc_distribution($fid);
print Dumper(\@gc_dist );

my @base_dist = EASIH::QC::db::fetch_base_distribution($fid);
print Dumper(\@base_dist );
