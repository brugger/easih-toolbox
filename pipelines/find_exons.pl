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

use lib '/software/lib/e62/ensembl/modules/';
use lib '/software/lib/bioperl/';

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 0;
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


use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;

my $species     = "human";
my $host        = 'mgpc17.medschl.cam.ac.uk';
my $user        = "easih_ro";
# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);

my $exit_count = 9;
# get adaptors
my $sa  = $reg->get_adaptor($species, 'core', 'slice');
my $ga  = $reg->get_adaptor($species, 'core', 'gene');
my $ta  = $reg->get_adaptor($species, 'core', 'transcript');

my $transcript_id  = 'NM_080680';
$transcript_id = shift;


my @genes = @{ $ga->fetch_all_by_external_name( $transcript_id )};

my @exons;

my $gene_name = "";
my $refseq_id = "";
my $strand    = 1;
foreach my $gene ( @genes ) {
  
  my $transcripts = $gene->get_all_Transcripts();
  while ( my $transcript = shift @{$transcripts} ) {
    my @db_entries = @{ $transcript->get_all_DBEntries('RefSeq_dna')};
    
    foreach my $db_entry (@db_entries) {

      $refseq_id = $db_entry->primary_id();
      $gene_name = $gene->external_name();
      if ($refseq_id eq $transcript_id ) {
	foreach my $exon ( @{ $transcript->get_all_translateable_Exons() } ) {


	  my $seq_region = $exon->slice->seq_region_name();
	  my $start      = $exon->start();
	  my $end        = $exon->end();
	  $strand     = $exon->strand();

	  unshift @exons, [$seq_region, $start, $end];
	}
      }
    }
  }
}

my $counter = 1;
$counter = int (@exons ) if ( $strand == -1 );
foreach my $exon ( @exons ) {
  print join("\t", @$exon, $gene_name ."_Exon". $counter, $transcript_id) . "\n";
  $counter++ if ( $strand == 1 );
  $counter-- if ( $strand == -1 );
}

# 
# 
# 
# Kim Brugger (09 Nov 2010)
sub usage {
  
  $0 =~ s/.*\///;

#  print "USAGE: $0 -b[am file] -i[indel vcf file] -s[np vcf file] -T<ranform, use if mapped against hg18> -B[ait file] -l[eeway, default 100 bp] -c[ount bases, need a -b as well]\n";
  print "USAGE: $0 -b[am file] -v[ariant vcf file] -T<ranform, use if mapped against hg18> -B[ait file] -l[eeway, default 100 bp] -c[ount bases, need a -b as well]\n";

#  print "\nor extrapolate the standard <bam, SNP vcf, indel vcf, output files> with the -Q <basefile name> option\n";
  print "\nor extrapolate the standard <bam, vcf, output files> with the -Q <basefile name> option\n";
  print "EXAMPLE: $0 -Q [base name] -T<ransform>\n";
  print "\n";

  
  print "USAGE: -o[output file]\n";  

  exit;
}



