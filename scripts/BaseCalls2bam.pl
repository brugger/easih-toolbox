#!/usr/bin/perl 
# 
# Transforms a illumina basecalls folder into 8 or 16 fq files depending on 
# is being a single or a paired ends run. 
# 
# Kim Brugger (10 Aug 2010), contact: kim.brugger@easih.ac.uk

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
use EASIH::DONE;
use EASIH::Sample;
use EASIH::Illumina::Summary;
use EASIH::Illumina::Sample_sheet;
use EASIH::Illumina::Config;
use EASIH::Parallel;
use EASIH::MD5;


my %opts;
#getopts("a:A:1:2:3:4:5:6:7:8:hs:Si:o:lhmb:L:", \%opts);
getopts("hB:L:bfSdo:D", \%opts);


# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub usage {
  $0 =~ s/.*\///;
  print "USAGE: $0 extracts data from the illumina BaseCalls directory.\n";
  print "USAGE: -h<elp>\n";
  print "USAGE: -i<nput dir, default BaseCalls>.\n";
  print "USAGE: -s<ample shee, normally picked up from cwd or BaseCalls>.\n";

  print "\nUSAGE: Advanced parameters\n";
  print "USAGE: -d<ata monger mode, oh dont you dare>\n";
  print "USAGE: -b<arcoded run>.\n";
  print "USAGE: -l<imited lanes, by default the whole slide is extracted>.\n";
  print "USAGE: -m<ismatches 1 bp error in barcodes>.\n";
  print "USAGE: -o<utput dir, default is /data/<project ID>/raw/, on a sample basis>.\n";
  print "USAGE: for barcoded lanes, please use a sample sheet.\n";
  print "USAGE: naming a (or all) lane with switches overrules all sample sheet checking.\n";

  exit -1;
}


my $basecall_folder   = $opts{'B'} || "./";
# I need to find the run_folder name, for the data mongering things.
my $run_folder   = id_run_folder();
my $intensities_folder = "$basecall_folder/../";
# The config file is checked to see if this was an indexed run anyway
my $run_config = EASIH::Illumina::Config::readin( $basecall_folder );
my $indexed_run = 0;
$indexed_run = 1 if ( $$run_config{ is_indexed });
my $reads = @{$$run_config{ reads}};
my $PE = 0;
$PE = 1 if ( ( $indexed_run && $reads >= 3) ||
	     (!$indexed_run && $reads == 2));

my @lanes = @{$$run_config{ 'lanes'}};
@lanes = verify_lane_names( \@lanes, split(',', $opts{'L'})) if ( $opts{'L'} ) ;
#print "Indexed: $indexed_run || PE: $PE || Lanes: @lanes || Reads: $reads\n";
my $sample_sheet = $opts{'s'} || find_sample_sheet( $basecall_folder );


my $mismatches    = $opts{'m'} || 0;
my $datamonger  = $opts{'d'} || 1;
my $outdir      = $opts{'o'} || "/data/";
my $parallel =  1;
$parallel = 0 if ($opts{S});

my $paired_data = 0;

my $fq_out  = $opts{ 'f' } || 1;
my $bam_out = $opts{ 'b' } || 0;

my $debug = $opts{ 'D' } || 0;
#$debug = 1;
EASIH::DONE::Connect('done_dev') if ($debug); 
$outdir = "/tmp/BCL2BAM/" if ($debug);



usage() if (! $sample_sheet || $opts{h});

my $rid = EASIH::DONE::add_run($run_folder, 'ILLUMINA');
print "RID :: $rid \n";
 

my %filenames;
my $sample_names = validate_names_and_open_outfiles( $sample_sheet, @lanes );

print Dumper($sample_names);


#$runfolder = "ILL_TEST_10" if ( $debug );

#my $res = EASIH::Illumina::Summary::readin_summaries($basecall_folder);
#EASIH::DONE::add_illumina_lane_stats_summary( $rid,  $res );

#exit;


fail("no sample sheet!\n", "MISSING_SAMPLESHEET") if ($datamonger && ! $sample_sheet || ($sample_sheet && ! -e $sample_sheet));



my (%fhs, %fids);

print Dumper( \%filenames );

#exit;

foreach my $lane_nr ( @lanes ) {
#for(my $lane = 4; $lane<=8; $lane++) {
  next if (!$$sample_names{ $lane_nr });
  EASIH::Parallel::job_push(\&analyse_lane, $lane_nr );
}


if ( $parallel ) {
  print EASIH::Parallel::run_parallel( );
}
else {
  EASIH::Parallel::run_serial( );
}

foreach my $lane_nr  (keys %filenames) {
  # Close all the file handles to the fq files.
  foreach my $fh ( keys %{$fhs{ $lane_nr }} ) {
    close ($fhs{ $lane_nr }{$fh}) || die "Could not close filehandle '$fh': $!\n";
  }
  # queue up the calculation of md5 sums
  foreach my $filename (keys %{$filenames{ $lane_nr }}) {
    EASIH::Parallel::job_push(\&EASIH::MD5::create_file, "$filenames{ $lane_nr }{$filename}");
  }
}


if ( $parallel ) {
  print EASIH::Parallel::run_parallel( );
}
else {
  EASIH::Parallel::run_serial( );
}

if ( $datamonger ) {
    my $res = EASIH::Illumina::Summary::readin_summaries($basecall_folder);
    EASIH::DONE::add_illumina_lane_stats_summary( $rid,  $res );
    EASIH::DONE::add_offloading_status($rid, 
				       "BASECALLS2FQ_DONE");
}



# 
# 
# 
# Kim Brugger (17 Jun 2011)
sub outfile {
  my ( $barcode, $lane_nr ) = @_;

#  return -1;

#  print "$barcode, $lane_nr \n";

  return $fhs{ $lane_nr }{ $barcode } if ( $fhs{ $lane_nr }{ $barcode } );

#  my ($basename, $read_nr) = $sample_name =~ /^(.*?)\.[12]/;

  my $sample_filename = $filenames{ $lane_nr }{ $barcode };

#  print "$lane_nr $barcode $sample_filename\n";
  
  if ( $sample_filename =~ /\.1\.fq\.gz/ ||  $sample_filename =~ /\.2\.fq\.gz/ ) {
    my $base_filename = $filenames{ $lane_nr }{ $barcode } || die "No filename for '$sample_filename' in lane '$lane_nr'\n";
    
    my $fh;
    open ($fh, "| gzip -c > $sample_filename") || fail( "Could not open '$base_filename': $!\n", "BASECALL2FQ_PATH_ERROR");
    $fhs{ $lane_nr }{ $barcode }  = $fh;
  }
  elsif ( $sample_filename =~ /\.bam/ ) {
    my $base_filename = $filenames{ $lane_nr }{ $barcode } || die "No entry in filenames for '$sample_filename' in lane '$lane_nr'\n";
    
    my $fh;
    open ($fh, "| samtools view -Sb - > $sample_filename") || fail( "Could not open '$base_filename': $!\n", "BASECALL2FQ_PATH_ERROR");
    $fhs{ $lane_nr }{ $barcode }  = $fh;
  }
  
  if ( $datamonger ) {
    
    my ($sample, $project) = EASIH::Sample::filename2sampleNproject($sample_filename);
    my $fid = EASIH::DONE::add_file($sample_filename, $sample, $project, $run_folder, 'ILLUMINA');
    $fids{ $lane_nr }{ $sample_filename }  = $fid;
  }
  
#  return $fh;
}

# 
# 
# 
# Kim Brugger (29 Oct 2012)
sub validate_names_and_open_outfiles {
  my ( $sample_sheet, @lanes ) = @_;

  if ( $outdir ) {
    if ( -e "$outdir" && ! -d "$outdir") {
      fail("$outdir is not a directory\n", "BASECALL2FQ_PATH_ERROR");
    }
    elsif ( -e "$outdir" && ! -w "$outdir") {
      fail("$outdir is not writeable\n", "BASECALL2FQ_PATH_ERROR");
    }
    if ( $outdir && ! -d $outdir ) {
      system "mkdir -p $outdir" || die "Could not create directory '$outdir': $!\n";
    }
  }

  my ($sample_names, $errors, $warnings, $project_name ) = EASIH::Illumina::Sample_sheet::readin( $sample_sheet );
#  print Dumper( $sample_names );

  # if there is only one barcode in a lane change this so the whole lane gets assigned to the sample in the lane.
  foreach my $lane_nr ( keys %$sample_names ) { 
    if ( keys %{$$sample_names{ $lane_nr }} == 1 ) {
      my ($barcode) = keys %{$$sample_names{ $lane_nr }};
      if ( $barcode ne "default" && $barcode ne "") {
	$$sample_names{ $lane_nr }{"default"} = $$sample_names{ $lane_nr }{$barcode};
	delete $$sample_names{ $lane_nr }{$barcode};
      }
    }
  }

#  die Dumper( $res );
  
#  print Dumper( $sample_names );
#  print Dumper( \%filenames );

  my %DONE_sample_sheet = EASIH::DONE::fetch_sample_sheet_hash( $rid );
#  print STDERR Dumper( \%DONE_sample_sheet );
  foreach my $lane_nr ( keys %DONE_sample_sheet ) {
    my %barcodes;
    map { s/\..*//; $barcodes{$_}++} keys %{$DONE_sample_sheet{ $lane_nr }};
    if (keys %barcodes == 1 ) {
      foreach my $barcode ( keys %{$DONE_sample_sheet{ $lane_nr }} ) {
	my ($post_fix) = $barcode =~ /\.(.*)/;
	$DONE_sample_sheet{ $lane_nr }{"default.$post_fix"} = $DONE_sample_sheet{ $lane_nr }{ $barcode };
	delete $DONE_sample_sheet{ $lane_nr }{ $barcode } if ($barcode ne "default.$post_fix");
      }
    }
  }
  
#  print STDERR Dumper( \%DONE_sample_sheet );
#  print Dumper( $sample_names );


#  fail( $errors, "MALFORMED_SAMPLESHEET" ) if ($errors);

  $errors = EASIH::Illumina::Sample_sheet::validate( $sample_names, @lanes ) ;
#  fail( $errors, "MALFORMED_SAMPLESHEET") if ( $errors );

  # assign filenames to each sample in each lane, as this script can/will
  # run in parallel this has to be done before we loose control.
  foreach my $lane_nr ( @lanes ) {
    next if ( ! $$sample_names{ $lane_nr });

    foreach my $bcode (keys %{$$sample_names{$lane_nr}}) {

      my $sample_name = $$sample_names{ $lane_nr }{ $bcode };

#      print STDERR "LANE $lane_nr $bcode $sample_name\n";
      
      my ($base_filename, $error);

      if ( $project_name && $project_name =~ /^CP\d+/) {
	($base_filename, $error) = EASIH::Sample::sample2outfilename_wo_project_dir_n_version( "$sample_name", "$outdir/CP/$project_name/");
	system "cp $sample_sheet /data/CP/$project_name/";
      }
      else {
	($base_filename, $error) = EASIH::Sample::sample2outfilename( "$sample_name", $outdir);
      }

#      print STDERR "LANE $lane_nr $bcode $sample_name $base_filename\n";

      if ( $fq_out ) { 
      
	$filenames{ $lane_nr }{ "$bcode.1" } = "$base_filename.1.fq.gz";

	$filenames{ $lane_nr }{ "$bcode.1" } = $DONE_sample_sheet{ $lane_nr }{ "$bcode.1" } 
	  if ($DONE_sample_sheet{ $lane_nr }{ "$bcode.1" });

	outfile("$bcode.1", $lane_nr);
	if ( $PE ) {
	  $filenames{ $lane_nr }{ "$bcode.2" }   = "$base_filename.2.fq.gz";
	  $filenames{ $lane_nr }{ "$bcode.2" } = $DONE_sample_sheet{ $lane_nr }{ "$bcode.2" } 
	    if ($DONE_sample_sheet{ $lane_nr }{ "$bcode.2" });
	  outfile("$bcode.2", $lane_nr) if ( $PE && $fq_out);
	}
      }
      
      if ( $bam_out ) {
	$filenames{ $lane_nr }{ "$bcode.bam" } = "$base_filename.bam";
	$filenames{ $lane_nr }{ "$bcode.bam" } = $DONE_sample_sheet{ $lane_nr }{ "$bcode.bam" } 
	  if ($DONE_sample_sheet{ $lane_nr }{ "$bcode.bam" });
	outfile("$bcode.bam", $lane_nr);
      }
    }
  }

#  die Dumper( \%filenames );
#  print Dumper ( \%filenames );
  
  return $sample_names;
}



# 
# 
# 
# Kim Brugger (06 Jan 2011)
sub verify_bcode {
  my ($bc1, @bc2s) = @_;

  return 'default' if ( $bc1 eq 'default');
  return 'default' if ( $bc2s[0] eq 'default');
  
  if (length($bc1) > length($bc2s[0])) {
    $bc1 = substr( $bc1, 0, length($bc2s[0]));
  }

#  exit;
  foreach my $bc2 ( @bc2s ) {

    if ( $mismatches == 0) {
      return $bc1 if ($bc1 eq $bc2);
      next;
    }
    else {
      my @seq1 = split('', $bc1);
      my @seq2 = split('', $bc2);
      
      my $diffs = 0;
      for( my $i = 0; $i < @seq1 && $i < @seq2; $i++) {
	$diffs++ if ( $seq1[$i] ne $seq2[$i]);
	last if ( $diffs > 1);
	
      }
      
      return $bc2 if ( $diffs <= 1 );
    }
  }

  return undef;
}




# 
# Need the run_folder for the datamongering. As the script can be
# called in every possible way this is a tad complicated
# 
# Kim Brugger (17 Jun 2011)
sub id_run_folder {

  my $dir = "";

  our $cwd      = `pwd`;
  chomp($cwd);


  # An absolute path was given to the BaseCalls dir
  if ( $basecall_folder && $basecall_folder =~ /^\// ) {
    $dir = $basecall_folder;
  }
  elsif( $basecall_folder ) {
    $dir = "$cwd//$basecall_folder";
  }
  else { 
    $dir = $cwd;
  }

  print "$dir $cwd\n";

  # the user might only have given the top run folder as input, so lets go on a guessing expedition
  $dir = "$dir/Data/" if ( -e "$dir/Data");
  $dir = "$dir/Intensities/" if ( -e "$dir/Intensities/");
  $dir = "$dir/BaseCalls/" if ( -e "$dir/BaseCalls");


  # remove double // in the name
  $dir =~ s/\/{2,}/\//g;
  # and any /./
  $dir =~ s/\/\.\//\//g;

  if ($dir !~ /Data\/Intensities\/BaseCalls/) {
    fail("Not a illumin runfolder/BaseCalls directory \n", "BASECALL2FQ_INPATH_ERROR");
  }
  

  if ( $dir  =~ /\.\./ ) {
    fail("Cannot handle input paths containing: ../\n", "BASECALL2FQ_INPATH_ERROR");
  }

  $basecall_folder = $dir;

  my @dirs = split( "/", $dir);
  return $dirs[3];
}


# 
# 
# 
# Kim Brugger (17 Jun 2011)
sub fail {
  my ($message, $status) = @_;

  
  if ( $datamonger ) {
    print STDERR "$message\n";
    EASIH::DONE::add_offloading_status($rid, 
				       "$status");
    EASIH::DONE::add_offloading_status($rid, 
				       "BASECALLS2FQ_FAILED");
    exit -1;
  }
  else {
    die $message;
  }
  
}



# 
# 
# 
# Kim Brugger (29 Oct 2012)
sub verify_lane_names {
  my ( $config_names, @param_names) = @_;

  my %config_names_hash;
  map { $config_names_hash{ $_ }++} @$config_names;

  map { fail( "$_ is not a valid lane name for this run\n", "INVALID_LANE_ID") if (! $config_names_hash{ $_ }) } @param_names;

  return @param_names;
}






# 
# find the runs sample sheet, as it can be named in a million ways, this is the solution for now.
# 
# Kim Brugger (29 Oct 2012)
sub find_sample_sheet {
  my ( $indir ) = @_;

  return "$indir/sample_sheet.csv" if ( -e "$indir/sample_sheet.csv");
  return "$indir/Sample_sheet.csv" if ( -e "$indir/Sample_sheet.csv");
  return "$indir/sample_Sheet.csv" if ( -e "$indir/sample_Sheet.csv");
  return "$indir/Sample_Sheet.csv" if ( -e "$indir/Sample_Sheet.csv");

  return "$indir/Samplesheet.csv" if ( -e "$indir/Samplesheet.csv");
  return "$indir/sampleSheet.csv" if ( -e "$indir/sampleSheet.csv");
  return "$indir/samplesheet.csv" if ( -e "$indir/samplesheet.csv");
  return "$indir/SampleSheet.csv" if ( -e "$indir/SampleSheet.csv");

  return "$indir/../../../sample_sheet.csv" if ( -e "$indir/../../../sample_sheet.csv");
  return "$indir/../../../Sample_sheet.csv" if ( -e "$indir/../../../Sample_sheet.csv");
  return "$indir/../../../sample_Sheet.csv" if ( -e "$indir/../../../sample_Sheet.csv");
  return "$indir/../../../Sample_Sheet.csv" if ( -e "$indir/../../../Sample_Sheet.csv");

  return "$indir/../../../Samplesheet.csv" if ( -e "$indir/Samplesheet.csv");
  return "$indir/../../../sampleSheet.csv" if ( -e "$indir/sampleSheet.csv");
  return "$indir/../../../samplesheet.csv" if ( -e "$indir/samplesheet.csv");
  return "$indir/../../../SampleSheet.csv" if ( -e "$indir/SampleSheet.csv");

  return undef;
}




# 
# 
# 
# Kim Brugger (30 Oct 2012)
sub analyse_lane {
  my ( $lane_nr ) = @_;

  my @bar_codes = sort keys %{$$sample_names{ $lane_nr }} if (ref $$sample_names{ $lane_nr } eq "HASH");
#  print "@bar_codes\n";

  my %stats;
  
  my $cmd = "illumina2bam -T Illumina2bam.jar I=$intensities_folder R=$basecall_folder/../../../ B=$basecall_folder L=$lane_nr O=/dev/stdout SC=EASIH  PF_FILTER=false | samtools view -h - | ";

  $cmd = "illumina2bam -T Illumina2bam.jar TILE_LIMIT=1 I=$intensities_folder R=$basecall_folder/../../../ B=$basecall_folder L=$lane_nr O=/dev/stdout SC=EASIH  PF_FILTER=false | samtools view -h - | " if ( $debug );

  print STDERR "$cmd\n";
#  exit;
  # The stats is done on a lane basis as we are running this bastard in parallel this is the way to do it. 
  # Way nicer than named fork-pipes etc.
  my %sample_stats;
  my $read_length = 0;

  my $prev_read_barcode = undef;
  my $strand = 1;

#  my $exit_counter = 5000;
  my ($lane_total, $lane_pass, $lane_bases, $lane_QV30) = (0,0,0,0);

  open ( my $i2b_pipe,  "$cmd") || die "Could not open illumina2bam pipe: $!\n";
  while (my $line = <$i2b_pipe>) {
    
    if ($line =~ /^\@/ ) {
      if ( $bam_out ) {
	foreach my $bcode ( keys %{$fhs{ $lane_nr }} ) {
	  if ($bcode =~ /\.bam/) {
	    my $fh = outfile($bcode, $lane_nr);
	    print $fh $line;
	  }
	}
      }
      next;
    }

    $lane_total++;
     
    my ($id, $flag, undef, undef, undef, undef, undef, undef, undef, $seq, $qual, $read_barcode) = split("\t", $line, 12);

    # read did not pass filter (PF == 0)
    next if ($flag & 0x0200 );
    $lane_pass++;

    $read_length = length( $seq ) if ( ! $read_length );

#    print "$read_barcode $strand\n";

    if ($flag & 0x0040 ) {
      $strand = 1;
    }
    if ($flag & 0x0080 ) {
      $strand = 2;
    }
 
    if ($indexed_run && $read_barcode =~ /.*?BC:Z:(.*?)\t.*/) {
      $read_barcode =~ s/.*?BC:Z:(.*?)\t.*/$1/;
      chomp( $read_barcode );
      $prev_read_barcode = $read_barcode;
    }
    else {
      $read_barcode = $prev_read_barcode;
    }

    if ( $indexed_run ) {
    
#      print "$read_barcode\n";
      $read_barcode = verify_bcode($read_barcode, @bar_codes) if ( $indexed_run );
      
      next if (! $read_barcode );
    }
      
    $read_barcode ||= 'default';

    my $sample_name = $$sample_names{ $lane_nr }{ $read_barcode } if ($$sample_names{ $lane_nr }{ $read_barcode });



    map{ $lane_QV30++ if (ord($_) - 33 > 30)} split("", $qual);

    $lane_bases += $read_length;
    
    if ( $fq_out ) {
      my $bfout = outfile("$read_barcode.$strand", $lane_nr);
      print $bfout "\@$id\n$seq\n+\n$qual\n";
      $sample_stats{"$read_barcode"}{"$strand"}++;
    }

    if ( $bam_out ) {
      my $bfout = outfile("$read_barcode.bam", $lane_nr );
      print $bfout $line;
      $sample_stats{"$read_barcode"}{"bam"}++;
    }

#    last if ( ! $exit_counter--);
  }


  if ( $datamonger ) {

    EASIH::DONE::add_illumina_lane_stats( $rid, 1, $lane_nr, $lane_total, $lane_pass, $lane_bases, $lane_QV30 );

    foreach my $barcode ( keys %sample_stats ) {
      foreach my $read ( keys %{$sample_stats{ $barcode }}) {
	next if ($read eq "bam");
	

	my $filename = $filenames{$lane_nr}{"$barcode.$read"};
	
	my $fid = $fids{ $lane_nr }{ $filename };
	

#print STDERR Dumper( \%fids );
#print STDERR Dumper( \%filenames );
#print STDERR Dumper( $sample_names );
#	exit;

	my $sample_name = $$sample_names{ $lane_nr }{ $barcode };
#	exit;
	
	print STDERR  "'$barcode.$read' --> $filename --> $fid\n";
	
	$barcode = "" if ( $barcode eq "default" );
	EASIH::DONE::add_illumina_sample_stats( $rid, $fid, $lane_nr, $read, $sample_name, $barcode);
	EASIH::DONE::update_file($fid , $sample_stats{$barcode}{$read}, $read_length);

	$barcode = "default" if ( $barcode eq "" );
	
	
      }
    }
  }

  
  
}



