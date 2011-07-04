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
my $DYNAMIC_LIB_PATHS = 1;
BEGIN {
  if ( $DYNAMIC_LIB_PATHS ) {
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
  else {
    use lib '/home/kb468/easih-toolbox/modules/';
  }
}

use EASIH;
use EASIH::QC;
use EASIH::QC::db;
use EASIH::Misc;
use EASIH::Sample;

my %opts;
getopts('b:c:f:q:hs:p:r', \%opts);

my $fastq_file    = $opts{f};
my $csfasta_file  = $opts{c};
my $qual_file     = $opts{q};
my $bam_file      = $opts{b};
my $sample_size   = $opts{s} || 20;
my $platform      = uc($opts{p}) || usage();
my $random_sample = $opts{r} || 0;
#sample limit MB
EASIH::QC::sample_size($sample_size);
EASIH::QC::random_sample(1) if ( $random_sample );


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
  print "USAGE: $0 -f[astq file] -c[sfasta file] -q[ual file] -b[am file] -p[latform, either SOLID or ILLUMINA or TORRENT] -s<ample size (in Mbases)> -r<andom sampling>\n";
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
	   "110628_MGILLUMINA2_00039_FC");


if ( $fastq_file || $csfasta_file || $qual_file ) {

  if ( $fastq_file ) {

    my $sname = `zcat -f $fastq_file | head -n1`;
    chomp $sname;
    $sname =~ s/\@(.*?_)(\d+).*/$1.sprintf("%05d", $2)/e;

    my $runfolder = (grep(/$sname/,  @orf))[0];
    $runfolder ||= $sname;

    my ($sample, $project) = EASIH::Sample::filename2sampleNproject($fastq_file);
    my $fid = EASIH::QC::db::add_file($fastq_file, $sample, $project, $runfolder, $platform);
    EASIH::QC::fid($fid);

    my ( $name, $total_reads, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC) = EASIH::QC::db::fetch_file_info( $fid );
    
    if ( ! $total_reads ) {
      my $fq_lines = `zcat -f $fastq_file | wc -l`;
      chomp $fq_lines;
      $fq_lines /= 4;
      print "$fq_lines in $fastq_file\n";
      EASIH::QC::db::update_file($fid, $fq_lines, undef, undef, undef, undef, undef);
    }


    $QC = EASIH::QC::fastQC( $fastq_file, 1 );

#    print Dumper( $QC );

    EASIH::QC::base_qual2db( $QC );
    EASIH::QC::base_dist2db( $QC );
    EASIH::QC::base_qual_dist2db( $QC );
    EASIH::QC::GC2db( $QC );
    EASIH::QC::duplicates2db( $QC );
    EASIH::QC::partial_adaptor2db( $QC );

    

  }

  if ( $csfasta_file ) {
    die "does not do csfasta files\n";
    $QC = EASIH::QC::csfastaQC( $csfasta_file);
    $base_name = $csfasta_file;
    $infile    = $csfasta_file;
    if ( $base_name =~ /\.\// || $base_name !~ /^\//) {
      $base_name =~ s/\.\///;
      my $cwd = `pwd`;
      chomp $cwd;
      $base_name = "$cwd/$base_name";
    }
    $base_name =~ s/(.*?)\..*/$1/;
    $base_name = "$cwd/$base_name" if ( $base_name !~ /\//);
    EASIH::QC::make_plots($QC, "$tmp_dir/$tmp_file", $base_name);
  }

  if ( $qual_file ) {
    die "does not do qual files\n";
    $QC = EASIH::QC::qualQC( $qual_file );


    $base_name = $qual_file;
    $infile    = $qual_file;
    if ( $base_name =~ /\.\// || $base_name !~ /^\//) {
      $base_name =~ s/\.\///;
      my $cwd = `pwd`;
      chomp $cwd;
      $base_name = "$cwd/$base_name";
    }
    $base_name =~ s/(.*?)\..*/$1/;
    $base_name = "$cwd/$base_name" if ( $base_name !~ /\//);

    EASIH::QC::make_plots($QC, "$tmp_dir/$tmp_file", $base_name);
  }


}
elsif ( $bam_file ) {
  my ($QC1, $QC2) = EASIH::QC::bamQC( $bam_file );
  die "does not do bam files\n";

  $base_name = $bam_file;
  $infile    = $bam_file;
  $base_name =~ s/(.*?)\..*/$1/;
  $base_name = "$cwd/$base_name" if ( $base_name !~ /\//);
  
  EASIH::QC::make_plots($QC1, "$tmp_dir/$tmp_file", $base_name);

  open (my $out, ">  $tmp_dir/$tmp_file.tex") || die "Could not open outfile: $!\n";
  print $out latex_header();
  print $out latex_summary($infile, uc($platform), $sample_size, $$QC1{reads} || $$QC1{quals}, $$QC1{Q30}, $$QC1{perc_dup}, $$QC1{mappability}, $$QC1{ACsplit}, $$QC{partial_adaptor} );
  print $out latex_QV("$tmp_dir/$tmp_file\_BaseQual.pdf", "$tmp_dir/$tmp_file\_QualHist.pdf");
  print $out latex_dups("$tmp_dir/$tmp_file\_DupHist.pdf", $$QC1{duplicates});
  print $out latex_pam("$tmp_dir/$tmp_file\_PAM.pdf", $$QC{partial_adaptor});
  print $out latex_GC("$tmp_dir/$tmp_file\_BaseDist.pdf", "$tmp_dir/$tmp_file\_GC.pdf");
  print $out latex_tail();
  make_pdf($tmp_dir, "$tmp_file.tex", "$base_name.1.pdf");
  system "rm -rf $tmp_dir/*";
  print "Report: $base_name.1.pdf\n";

  EASIH::QC::make_plots($QC2, "$tmp_dir/$tmp_file", $base_name);
  open ( $out, ">  $tmp_dir/$tmp_file.tex") || die "Could not open outfile: $!\n";
  print $out latex_header();
  print $out latex_summary($infile, uc($platform), $sample_size, $$QC2{reads} || $$QC2{quals}, $$QC2{Q30}, $$QC2{perc_dup}, $$QC2{mappability}, $$QC2{ACsplit}, $$QC{partial_adaptor} );
  print $out latex_QV("$tmp_dir/$tmp_file\_BaseQual.pdf", "$tmp_dir/$tmp_file\_QualHist.pdf");
  print $out latex_dups("$tmp_dir/$tmp_file\_DupHist.pdf", $$QC2{duplicates});
  print $out latex_pam("$tmp_dir/$tmp_file\_PAM.pdf", $$QC{partial_adaptor});
  print $out latex_GC("$tmp_dir/$tmp_file\_BaseDist.pdf", "$tmp_dir/$tmp_file\_GC.pdf");
  print $out latex_tail();
  make_pdf($tmp_dir, "$tmp_file.tex", "$base_name.2.pdf");
  system "rm -rf $tmp_dir";
  print "Report: $base_name.2.pdf\n";
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
