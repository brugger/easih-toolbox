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

my $debug = 0;

my %opts;
getopts('b:c:q:hs:p:r1:2:', \%opts);

my $fastq_file    = $opts{f} || $opts{1};
my $fastq_file2   = $opts{2};
die "-1 $fastq_file and -2 $fastq_file2 points to the same file\n" if ( $fastq_file && $fastq_file2 && 
									$fastq_file eq $fastq_file2 );
my $platform      = uc($opts{p}) || usage();
my $random_sample = $opts{r} || 0;
#sample limit MB
EASIH::QC::sample_size(-1);
EASIH::QC::sample_reads(100000);
EASIH::QC::sample_reads(50000) if ( $fastq_file2 );

EASIH::QC::sample_reads(1000) if ( $debug );
EASIH::QC::sample_reads(500) if ( $debug && $fastq_file2 );


EASIH::QC::random_sample(1) if ( $random_sample );
EASIH::QC::do_mappings(1);
EASIH::DONE::Connect('done_dev') if ($debug); 

my ( $tmp_dir, $tmp_file) = EASIH::Misc::tmp_dir_file();

usage() if ( !$fastq_file || ! $platform);
my $cwd      = `pwd`;
chomp($cwd);
my ($QC, $base_name, $infile);

$platform = "TORRENT" if ( $platform eq "ION");

# 
# 
# Kim Brugger (03 Feb 2011)
sub usage {
  $0 =~ s/.*\///;
  print "USAGE: $0 -1[fastq file] -2[fastq file] -p[latform, eitherILLUMINA or TORRENT]  -r<andom sampling>\n";
  exit -1;
}

if ( $fastq_file ) {

  die "$opts{1} and $opts{2} points to the  same file\n" if ( $fastq_file && $opts{2} && 
							      $fastq_file eq $opts{2} );
  
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
    exit if ( $fq_lines == 0);
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
  
  
  print Dumper( $QC1);
  
  
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
