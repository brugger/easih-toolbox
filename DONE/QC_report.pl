#!/usr/bin/perl 
# 
# Database version of the QC report script
# 
# 
# Kim Brugger (30 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

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

use EASIH;
use EASIH::QC;
use EASIH::DONE;
use EASIH::Misc;
use EASIH::Sample;

my %opts;
getopts('b:c:f:q:hs:p:r1:2:', \%opts);

my $fastq_file    = $opts{f} || $opts{1};
my $csfasta_file  = $opts{c};
my $qual_file     = $opts{q};
my $bam_file      = $opts{b};
my $sample_size   = $opts{s} || 20;
my $platform      = uc($opts{p}) || "ILLUMINA" || usage();
my $random_sample = $opts{r} || 0;
#sample limit MB
EASIH::QC::sample_size($sample_size);
EASIH::QC::random_sample(1) if ( $random_sample );
EASIH::QC::do_mappings(1);

my ( $tmp_dir, $tmp_file) = EASIH::Misc::tmp_dir_file();

usage() if ( ($fastq_file   && $csfasta_file) || 
	     ($fastq_file   && $qual_file) || 
	     ($fastq_file   && $bam_file) || 
#	     ($csfasta_file && $qual_file) || 
	     ($csfasta_file && $bam_file) || 
	     ($qual_file    && $bam_file) || 
	     ! $platform);


my $cwd      = `pwd`;
chomp($cwd);
my ($QC, $base_name, $infile);

$platform = "TORRENT" if ( $platform eq "ION");

my %platforms = ( "ILLUMINA" => "Illumina GA2X",
		  "SOLID"    => "SOLiD v4",
		  "TORRENT"  => "Ion Torrent");
# 
# 
# Kim Brugger (03 Feb 2011)
sub usage {
  $0 =~ s/.*\///;
  print "USAGE: $0 -1[fastq file] -2[fastq file] -p[latform, either SOLID or ILLUMINA or TORRENT] -s<ample size (in Mbases)> -r<andom sampling>\n";
  exit -1;
}

my @orf = ("101021_SOLEXA2_00011_FC",
	   "101021_SOLEXA2_00012_FC",
	   "101021_SOLEXA2_00013_FC",
	   "101022_SOLEXA2_00014_FC",
	   "101025_SOLEXA2_00015_FC",
	   "101108_SOLEXA2_00016_FC",
	   "101118_SOLEXA2_00017_FC_GW_FC4_181110",
	   "101126_SOLEXA2_00018_FC",
	   "101126_SOLEXA2_00019_FC_GW_FC5_261110",
	   "101207_SOLEXA2_00020_FC_PainFc6",
	   "101208_SOLEXA2_A06_FC6",
	   "101220_SOLEXA2_00022_FC_A23_A15_FC1_201210",
	   "110128_MGILLUMINA2_00024_FC",
	   "110405_MGILLUMINA2_00025_FC",
	   "110407_MGILLUMINA2_00026_FC_A06_FC7",
	   "110408_MGILLUMINA2_00027_FC_A06",
	   "110504_MGILLUMINA2_00028_FC_A35_A30",
	   "110505_MGILLUMINA2_00029_FC",
	   "110513_MGILLUMINA2_00030_FC_A04MH",
	   "110524_MGILLUMINA2_00031_FC_A33_FC2",
	   "110609_MGILLUMINA2_00032_FC",
	   "110609_MGILLUMINA2_00033_FC",
	   "110620_MGILLUMINA2_00034_FC",
	   "110624_MGILLUMINA2_00035_FC",
	   "110627_MGILLUMINA2_00036_FC",
	   "110627_MGILLUMINA2_00037_FC",
	   "110627_MGILLUMINA2_00038_FC",
	   "110628_MGILLUMINA2_00039_FC",
	   "110704_MGILLUMINA2_00040_FC",
	   "110715_MGILLUMINA2_00041_FC");


if ( $opts{1} || $fastq_file || $csfasta_file || $qual_file ) {

  if ( $fastq_file ) {


    die "$opts{1} and $opts{2} points to the  same file\n" if ( $fastq_file && $opts{2} && 
								$fastq_file eq $opts{2} );
    
    EASIH::QC::sample_reads( 1000000 );
    EASIH::QC::random_sample(1);
    my $fastq_file2 = $opts{2};

    my $sname = `zcat -f $fastq_file | head -n1`;
    chomp $sname;
    $sname =~ s/\@(.*?_)(\d+).*/$1.sprintf("%05d", $2)/e;

    my $runfolder = (grep(/$sname/,  @orf))[0];
    $runfolder ||= $sname;

    my ($sample, $project) = EASIH::Sample::filename2sampleNproject($fastq_file);
    my $fid1 = EASIH::DONE::add_file($fastq_file,  $sample, $project, $runfolder, $platform);
    my $fid2 = EASIH::DONE::add_file($fastq_file2, $sample, $project, $runfolder, $platform) if ( $fastq_file2);

    my ( $name, $total_reads, $read_length, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC) = EASIH::DONE::fetch_file_info( $fid1 );
    
    if ( ! $total_reads ) {
      my $fq_lines = `zcat -f $fastq_file | wc -l`;
      chomp $fq_lines;
      $fq_lines /= 4;
      print "$fq_lines in $fastq_file\n";
      EASIH::DONE::update_file($fid1, $fq_lines, undef, undef, undef, undef, undef, undef);
    }

    if (! $read_length ) {
      my $read_length = find_read_length( $fastq_file );
      EASIH::DONE::update_file($fid1, undef, $read_length, undef, undef, undef, undef, undef);
    }


    if ( $fid2 ) {
      my ( $name, $total_reads, $read_length, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC) = EASIH::DONE::fetch_file_info( $fid2 );
    
      if ( ! $total_reads ) {
	my $fq_lines = `zcat -f $fastq_file2 | wc -l`;
	chomp $fq_lines;
	$fq_lines /= 4;
	print "$fq_lines in $fastq_file2\n";
	EASIH::DONE::update_file($fid2, $fq_lines, undef, undef, undef, undef, undef, undef);
      }

      if (! $read_length ) {
	my $read_length = find_read_length( $fastq_file2 );
	EASIH::DONE::update_file($fid2, undef, $read_length, undef, undef, undef, undef, undef);
      }
    }

    my ( $QC1, $QC2 ) = EASIH::QC::fastQC( $fastq_file, $fastq_file2 );

    
    if ( $fid2 ) { 
      EASIH::QC::fid( $fid2 );
      EASIH::QC::base_qual2db( $QC2 );
      EASIH::QC::base_dist2db( $QC2 );
      EASIH::QC::base_qual_dist2db( $QC2 );
      EASIH::QC::GC2db( $QC2 );
      EASIH::QC::duplicates2db( $QC2 );
      EASIH::QC::partial_adaptor2db( $QC2 );
    }

    EASIH::QC::fid( $fid1 );
    EASIH::QC::base_qual2db( $QC1 );
    EASIH::QC::base_dist2db( $QC1 );
    EASIH::QC::base_qual_dist2db( $QC1 );
    EASIH::QC::GC2db( $QC1 );
    EASIH::QC::duplicates2db( $QC1 );
    EASIH::QC::partial_adaptor2db( $QC1 );
    
    EASIH::QC::mappings2db( $QC1, $fid2 );

  }
}



# 
# 
# 
# Kim Brugger (22 Jul 2011)
sub find_read_length {
  my ($name) = @_;
  
  my $f1;
  if ( $name =~ /gz/) {
    open ( $f1, "gunzip -c $name | ") || die "Could not open '$name': $!\n";
  }
  else {
    open ( $f1, "$name") || die "Could not open '$name': $!\n";
  }
  
  my $id   = <$f1>;
  my $seq  = <$f1>;
  my $str  = <$f1>;
  my $qual = <$f1>;

  close ($f1);
  chomp($seq);
  
  my $read_length = length($seq);

  return $read_length;
}




__END__

# 
# 
# 
# Kim Brugger (18 Aug 2010)
sub usage {
  
  $0 =~ s/.*\///;
  print "There are three running modes, one where the seq is sampled and aligned and the other from a pre aligned file, finally just qv and base dist.\n";
  print "USAGE RAW: $0 -1 [first fq file] -2 [second fq] -R[eference bwa formatted] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  print "USAGE RAW: $0 -b[am file] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  print "USAGE RAW: $0 -1 [first fq file] -2 [second fq] -n[o mapping] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  exit;

}
