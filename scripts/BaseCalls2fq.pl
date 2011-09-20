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

#my $debug = 0;
my $debug = 1;

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
use EASIH::Illumina::Sample_sheet;
use EASIH::Parallel;

my %opts;
getopts("a:A:1:2:3:4:5:6:7:8:hs:Si:o:lhnbd", \%opts);

my $limited_lanes = $opts{'l'};
my $no_mismatches = $opts{n};

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
  print "USAGE: -1,  -2, .., -8 < the lanes are assigned this name> Overrules the sample_sheet.csv\n";
  print "USAGE: -a< id, all lanes that are not named specifically gets this id>.\n";
  print "USAGE: -A< id, all lanes that are not named specifically gets this id with a _counter to separate them>.\n";
  print "USAGE: -d<ata monger mode, oh dont you dare>\n";
  print "USAGE: -b<arcoded run>.\n";
  print "USAGE: -l<imited lanes, by default the whole slide is extracted>.\n";
  print "USAGE: -n<o mismatches in barcodes, normal is 1 error>.\n";
  print "USAGE: -o<utput dir, default is /data/<project ID>/raw/, on a sample basis>.\n";
  print "USAGE: for barcoded lanes, please use a sample sheet.\n";
  print "USAGE: naming a (or all) lane with switches overrules all sample sheet checking.\n";

  exit -1;
}


my $indir       = $opts{'i'} || "./";
# I need to find the run_folder name, for the data mongering things.
my $runfolder   = id_run_folder();
my $datamonger  = $opts{'d'} || 0;
my $outdir      = $opts{'o'};
my $parallel =  1;
$parallel = 0 if ($opts{S});


my $paired_data = 0;

my $indexed_run = $opts{b} || 0;
my $sample_sheet = $opts{'s'};
$sample_sheet = "$indir/sample_sheet.csv" if (!$sample_sheet && -e "$indir/sample_sheet.csv");
$sample_sheet = "$indir/Sample_sheet.csv" if (!$sample_sheet && -e "$indir/Sample_sheet.csv");
$sample_sheet = "$indir/sample_Sheet.csv" if (!$sample_sheet && -e "$indir/sample_Sheet.csv");
$sample_sheet = "$indir/Sample_Sheet.csv" if (!$sample_sheet && -e "$indir/Sample_Sheet.csv");
if (!$sample_sheet && -e "BaseCalls/sample_sheet.csv") {
  $indir = "BaseCalls";
  $sample_sheet = "$indir/sample_sheet.csv";
}
if (!$sample_sheet && -e "sample_sheet.csv") {
  $indir = "./";
  $sample_sheet = "sample_sheet.csv";
}
fail("no sample sheet!\n", "MISSING_SAMPLESHEET") if ($datamonger && ! $sample_sheet || ($sample_sheet && ! -e $sample_sheet));


usage() if (! $sample_sheet &&  ! $opts{a} && ! $opts{A} && ! $opts{1} && ! $opts{2} && ! $opts{3} && ! $opts{4} && 
	   ! $opts{5} && ! $opts{6} && ! $opts{7} && ! $opts{8}  ||  $opts{h});

my ($sample_names, $removed_samples)  = readin_sample_sheet( $sample_sheet) if ($sample_sheet);

$$sample_names{1} = $opts{'1'} if ($opts{'1'});
$$sample_names{2} = $opts{'2'} if ($opts{'2'});
$$sample_names{3} = $opts{'3'} if ($opts{'3'});
$$sample_names{4} = $opts{'4'} if ($opts{'4'});
$$sample_names{5} = $opts{'5'} if ($opts{'5'});
$$sample_names{6} = $opts{'6'} if ($opts{'6'});
$$sample_names{7} = $opts{'7'} if ($opts{'7'});
$$sample_names{8} = $opts{'8'} if ($opts{'8'});

if ($opts{a}) {
  $$sample_names{1} = $opts{a} if ( !$opts{'1'} );
  $$sample_names{2} = $opts{a} if ( !$opts{'2'} );
  $$sample_names{3} = $opts{a} if ( !$opts{'3'} );
  $$sample_names{4} = $opts{a} if ( !$opts{'4'} );
  $$sample_names{5} = $opts{a} if ( !$opts{'5'} );
  $$sample_names{6} = $opts{a} if ( !$opts{'6'} );
  $$sample_names{7} = $opts{a} if ( !$opts{'7'} );
  $$sample_names{8} = $opts{a} if ( !$opts{'8'} );
}

if ($opts{A}) {
  my $counter = 1;
  $$sample_names{1} = "$opts{A}_".$counter++ if ( !$opts{'1'} );
  $$sample_names{2} = "$opts{A}_".$counter++ if ( !$opts{'2'} );
  $$sample_names{3} = "$opts{A}_".$counter++ if ( !$opts{'3'} );
  $$sample_names{4} = "$opts{A}_".$counter++ if ( !$opts{'4'} );
  $$sample_names{5} = "$opts{A}_".$counter++ if ( !$opts{'5'} );
  $$sample_names{6} = "$opts{A}_".$counter++ if ( !$opts{'6'} );
  $$sample_names{7} = "$opts{A}_".$counter++ if ( !$opts{'7'} );
  $$sample_names{8} = "$opts{A}_".$counter++ if ( !$opts{'8'} );
}

$sample_names = validate_lane_names( $sample_names);

my (%fhs, %fids);

my $rid = EASIH::DONE::add_run($runfolder, 'ILLUMINA') if ($datamonger);
my %reads_pr_sample;

#for(my $lane = 1; $lane<=8; $lane++) {
for(my $lane = 3; $lane<=8; $lane++) {

  next if (!$$sample_names{ $lane });

  EASIH::Parallel::job_push(\&analyse_lane, $lane);
#  analyse_lane($lane)
}


if ( $parallel ) {
  print EASIH::Parallel::run_parallel( );
}
else {
  EASIH::Parallel::run_serial( );
}


if ( $datamonger ) {

  # added the summed numbers to the file table, done as an update 
  foreach my $fid ( keys %reads_pr_sample ) {
    EASIH::DONE::update_file($fid, $reads_pr_sample{$fid});
  }

  EASIH::DONE::add_offloading_status($runfolder, 
				     "ILLUMINA", 
				     "BASECALLS2FQ_DONE");
}



# 
# 
# 
# Kim Brugger (17 Jun 2011)
sub open_outfile {
  my ( $sample_name, $lane_nr ) = @_;

  return $fhs{ $lane_nr }{ $sample_name } if ( $fhs{ $lane_nr }{ $sample_name } );

#  my ($basename, $read_nr) = $sample_name =~ /^(.*?)\.[12]/;

  use EASIH::Sample;
  my ($base_filename, $error) = EASIH::Sample::sample2outfilename( $sample_name, $outdir);

  # simplifies ////// to /
  $base_filename =~ s/\/{2,}/\//;
  $base_filename .= ".fq.gz";


  my $fh;
  open ($fh, "| gzip -c > $base_filename") || fail( "Could not open '$base_filename': $!\n", "BASECALL2FQ_PATH_ERROR");
  $fhs{ $lane_nr }{ $sample_name }  = $fh;
  
  if ( $datamonger ) {
    
    my ($sample, $project) = EASIH::Sample::filename2sampleNproject($base_filename);
    my $fid = EASIH::DONE::add_file($base_filename, $sample, $project, $runfolder, 'ILLUMINA');
    $fids{ $sample_name }  = $fid;
  }
  
  return $fh;
}


# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub analyse_lane {
  my ( $lane_nr ) = @_;

  # There is a need to go from sample name to barcode later on.
  my %bcode2sample;
  if ( ref($$sample_names{ $lane_nr }) eq 'HASH' && ! $$sample_names{ $lane_nr }{'default'}) {
    foreach my $bcode  ( keys %{$$sample_names{ $lane_nr }} ) {
      $bcode2sample{ $$sample_names{ $lane_nr }{ $bcode }} = $bcode;
    }    
  }
  
  # Get all the files for the lane.
  my @files = glob("$indir/s_$lane_nr\_1_*_qseq.txt");

  my ($lane_total, $lane_pass) = (0, 0);
  my %sample_stats;

  foreach my $file (@files) {
    my $tile_stats = analyse_tile( $file, $lane_nr );
    $lane_total += $$tile_stats{ 'lane_total' } || 0;
    $lane_pass  += $$tile_stats{ 'lane_pass'  } || 0;

    foreach my $k ( keys %$tile_stats ) {
      next if ( $k =~ /^lane_/ );
      $sample_stats{ $k } += $$tile_stats{ $k };
    }

    last if ($debug);
  }

  printf("lane $lane_nr\t\t$lane_total\t$lane_pass (%.2f %%)\t%.2f avg clusters per tile\n", $lane_pass*100/$lane_total, $lane_pass/120);

  foreach my $k ( keys %sample_stats ) {
    my $perc = sprintf("%.2f", $sample_stats{$k}*100/$lane_pass);
    my $barcode = $bcode2sample{$k} || "";
    printf("lane $lane_nr\t$k\t$barcode\t$sample_stats{$k}\t$perc %%\n");
  }

  if($datamonger) {
#    EASIH::DONE::add_illumina_lane_stats( $rid, $fids{"$sample.1"}, $lane_nr, 1, $sample, $lane_total, $lane_pass );

#     foreach my $k ( keys %{$multiplex_stats{pass}} ) {
#       my $sample_name = $$sample_names{ $lane_nr }{ $k };
#       my $perc = sprintf("%.2f", $multiplex_stats{pass}{$k}*100/$multiplex_stats{total});
#       printf("lane $lane_nr\t$sample_name\t$k\t$multiplex_stats{pass}{$k}\t$perc %%\n");
      
#       EASIH::DONE::add_illumina_multiplex_stats( $rid, $fids{"$sample_name.1"}, $lane_nr, 1, $sample_name, $k, $multiplex_stats{pass}{$k}, $multiplex_stats{pass}{$k} + ($multiplex_stats{fail}{$k} || 0), $perc);
      
#       $reads_pr_sample{$fids{"$sample_name.1"}} += $multiplex_stats{pass}{$k} if ( $datamonger);
#       $reads_pr_sample{$fids{"$sample_name.2"}} += $multiplex_stats{pass}{$k} if ( $datamonger);
#     }
  }


#  return ($in1, $out1, $in2, $out2);
}

# 
# 
# 
# Kim Brugger (15 Sep 2011)
sub analyse_tile {
  my ($file1, $lane_nr ) = @_;

  my ($file2, $file3) = ($file1, $file1);
  
  $file2 =~ s/(s_\d)_1_/$1_2_/;
  $file3 =~ s/(s_\d)_1_/$1_3_/;

  my $tile2seq = "/home/kb468/easih-toolbox/C/tile2seq/tile2seq";

  my ($fh1, $fh2, $fh3);
  open ($fh1, "$tile2seq $file1 |") || die "Could not open '$file1': $!\n" if ( -e $file1 );
  open ($fh2, "$tile2seq $file2 |") || die "Could not open '$file2': $!\n" if ( -e $file2 );
  open ($fh3, "$tile2seq $file3 |") || die "Could not open '$file3': $!\n" if ( -e $file3 );

  # this is a multiplexed run, so flip the file handle for file 2 & 3.
  # I am so going to regret this later on, I am sure of it...
  ($fh2, $fh3) = ($fh3, $fh2) if ( $fh3 );

  
  my %stats;
    
  if ( $fh3 ) {

    my @bar_codes = sort keys %{$$sample_names{ $lane_nr }} if (ref $$sample_names{ $lane_nr } eq "HASH");
    
    while (my $line3 = <$fh3>) {
      chomp ( $line3 );
      $stats{lane_total}++;
      
      my $line1 = <$fh1>;
      my $line2 = <$fh2> if ($fh2);
      my ($name3, $bases3, $qual3, $pf3) = split /\t/, $line3;
      next if ( ! $pf3 );

      chomp ( $line1 );
      chomp ( $line2 );

      my $barcode = verify_bcode($bases3, @bar_codes) if (ref $$sample_names{ $lane_nr } eq "HASH");
      my ($name1, $bases1, $qual1, $pf1) = split /\t/, $line1;
      my ($name2, $bases2, $qual2, $pf2) = split /\t/, $line2 if ($fh2);
      next if (! $barcode && ref $$sample_names{ $lane_nr } eq "HASH");

      my $sample_name = $$sample_names{ $lane_nr };
      $sample_name = $$sample_names{ $lane_nr }{ $barcode } if ( $barcode && ref $$sample_names{ $lane_nr } eq "HASH");



      my $bfout = open_outfile("$sample_name.1", $lane_nr );
      print $bfout "$name1\n$bases1\n+\n$qual1\n";

      if ($fh2) {
	my $bfout = open_outfile("$sample_name.2", $lane_nr );
	print $bfout "$name2\n$bases2\n+\n$qual2\n";
      }

      $stats{ $sample_name }++;
      $stats{lane_pass}++;
    }
  }
  else {
    my ($input_file, $fout, $demultiplexing, $ebcs) = @_;

    while (my $line1 = <$fh1>) {
      $stats{lane_total}++;
      
      my $line2 = <$fh2> if ($fh2);
      chomp ( $line1 );
      chomp ( $line2 );
      
      my ($name1, $bases1, $qual1, $pf1) = split /\t/, $line1;
      my ($name2, $bases2, $qual2, $pf2) = split /\t/, $line2 if ($fh2);

      next if ( ! $pf1 );

      my $sample_name = $$sample_names{$lane_nr};

      my $fout = open_outfile("$sample_name.1", $lane_nr );
      print $fout "$name1\n$bases1\n+\n$qual1\n";

      if ($fh2) {
	my $fout = open_outfile("$sample_name.2", $lane_nr );
	print $fout "$name2\n$bases2\n+\n$qual2\n";
      }

      $stats{ $sample_name }++;


      $stats{lane_pass}++;
    }
  }

  return (\%stats);
}

# 
# 
# 
# Kim Brugger (09 Jun 2011)
sub readin_sample_sheet {
  my ( $sample_sheet) = @_;
  
  my ($res, $errors ) = EASIH::Illumina::Sample_sheet::readin( $sample_sheet );

  fail( $errors, "MALFORMED_SAMPLESHEET" ) if ($errors);

  $indexed_run = EASIH::Illumina::Sample_sheet::indexed_run( $res );  


  ($res, my $removed_samples) = EASIH::Illumina::Sample_sheet::remove_easih_barcodes( $res );

#  die Dumper( $removed_samples );

  return ($res, $removed_samples);

  return %$res;
}


# 
# 
# 
# Kim Brugger (10 Aug 2010)
sub validate_lane_names {
  my ($sample_names) = @_;

  if ($sample_sheet ) {
    my $errors = EASIH::Illumina::Sample_sheet::validate( $sample_names, $limited_lanes ) ;
    fail( $errors, "MALFORMED_SAMPLESHEET") if ( $errors );
  }


  my %basenames;
  for ( my $lane =1; $lane <=8;$lane++) {
    
    next if ( ! $$sample_names{ $lane });

    if (ref ($$sample_names{$lane}) eq "HASH") {
      foreach my $bcode (keys %{$$sample_names{$lane}}) {
	$basenames{ $$sample_names{$lane}{$bcode}} = -1;
      }
    }
    else {
      $basenames{ $$sample_names{$lane} } = -1;
    }
  }
  


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


  return $sample_names;
}



# 
# 
# 
# Kim Brugger (06 Jan 2011)
sub verify_bcode {
  my ($bc1, @bc2s) = @_;

  foreach my $bc2 ( @bc2s ) {

    if ( $no_mismatches) {
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
# Need the rulfolder for the datamongering. As the script can be
# called in every possible way this is a tad complicated
# 
# Kim Brugger (17 Jun 2011)
sub id_run_folder {

  my $dir = "";

  our $cwd      = `pwd`;
  chomp($cwd);


  # An absolute path was given to the BaseCalls dir
  if ( $indir && $indir =~ /^\// ) {
    $dir = $indir;
  }
  elsif( $indir ) {
    $dir = "$cwd//$indir";
  }
  else { 
    $dir = $cwd;
  }

  # remove double // in the name
  $dir =~ s/\/{2,}/\//g;

  if ( $dir  =~ /\.\./ ) {
    fail("Cannot handle input paths containing: ../\n", "BASECALL2FQ_INPATH_ERROR");
  }

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
    EASIH::DONE::add_offloading_status($runfolder, 
				       "ILLUMINA", 
				       "$status");
    EASIH::DONE::add_offloading_status($runfolder, 
				       "ILLUMINA", 
				       "BASECALLS2FQ_FAILED");
    exit -1;
  }
  else {
    die $message;
  }
  
}

