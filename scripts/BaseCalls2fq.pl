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

my $debug = 0;
#my $debug = 1;

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

EASIH::Barcodes::barcode_set( 'ill09' );
EASIH::Barcodes::strict_tags();
EASIH::Barcodes::error_correct_barcodes(0);

my %opts;
getopts("a:A:1:2:3:4:5:6:7:8:hs:i:o:lhnbd", \%opts);

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

my %sample_names = readin_sample_sheet( $sample_sheet) if ($sample_sheet);


$sample_names{1} = $opts{'1'} if ($opts{'1'});
$sample_names{2} = $opts{'2'} if ($opts{'2'});
$sample_names{3} = $opts{'3'} if ($opts{'3'});
$sample_names{4} = $opts{'4'} if ($opts{'4'});
$sample_names{5} = $opts{'5'} if ($opts{'5'});
$sample_names{6} = $opts{'6'} if ($opts{'6'});
$sample_names{7} = $opts{'7'} if ($opts{'7'});
$sample_names{8} = $opts{'8'} if ($opts{'8'});

if ($opts{a}) {
  $sample_names{1} = $opts{a} if ( !$opts{'1'} );
  $sample_names{2} = $opts{a} if ( !$opts{'2'} );
  $sample_names{3} = $opts{a} if ( !$opts{'3'} );
  $sample_names{4} = $opts{a} if ( !$opts{'4'} );
  $sample_names{5} = $opts{a} if ( !$opts{'5'} );
  $sample_names{6} = $opts{a} if ( !$opts{'6'} );
  $sample_names{7} = $opts{a} if ( !$opts{'7'} );
  $sample_names{8} = $opts{a} if ( !$opts{'8'} );
}

if ($opts{A}) {
  my $counter = 1;
  $sample_names{1} = "$opts{A}_".$counter++ if ( !$opts{'1'} );
  $sample_names{2} = "$opts{A}_".$counter++ if ( !$opts{'2'} );
  $sample_names{3} = "$opts{A}_".$counter++ if ( !$opts{'3'} );
  $sample_names{4} = "$opts{A}_".$counter++ if ( !$opts{'4'} );
  $sample_names{5} = "$opts{A}_".$counter++ if ( !$opts{'5'} );
  $sample_names{6} = "$opts{A}_".$counter++ if ( !$opts{'6'} );
  $sample_names{7} = "$opts{A}_".$counter++ if ( !$opts{'7'} );
  $sample_names{8} = "$opts{A}_".$counter++ if ( !$opts{'8'} );
}


%sample_names = validate_lane_names(%sample_names);



my (%fhs, %fids);

my $rid = EASIH::DONE::add_run($runfolder, 'ILLUMINA') if ($datamonger);
my %reads_pr_sample;

for(my $lane = 1; $lane<=8; $lane++) {

  next if (!$sample_names{ $lane });
 
  # this lane is barcoded...
  if (ref ($sample_names{$lane}) eq "HASH" && ! $sample_names{$lane}{'default'}) {
    analyse_barcoded_lane($lane);
  }
  else {
    analyse_lane($lane)
  }
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
  my ($filename) = @_;

  

  # simplifies ////// to /
  $filename =~ s/\/{2,}/\//;

  my $fh;
#  open ($fh, "| gzip -c > $filename") || fail( "Could not open '$filename': $!\n", "BASECALL2FQ_PATH_ERROR");
  $filename =~ s/.gz//;
  open ($fh, " > $filename") || fail( "Could not open '$filename': $!\n", "BASECALL2FQ_PATH_ERROR");


  if ( $datamonger ) {
    
    my ($sample, $project) = EASIH::Sample::filename2sampleNproject($filename);
    my $fid = EASIH::DONE::add_file($filename, $sample, $project, $runfolder, 'ILLUMINA');
    $filename =~ s/.*\///;
    $filename =~ s/.fq.gz//;
    $filename =~ s/_\d+//;
    $fids{ $filename }  = $fid;
  }
  
  return $fh;
}


# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub analyse_lane {
  my ( $lane_nr ) = @_;

  my $sample_name = $sample_names{ $lane_nr }{'default'} if ( ref($sample_names{ $lane_nr }) eq 'HASH');
  $sample_name ||= $sample_names{ $lane_nr } ;
  my $basename = $sample_names{ $sample_name };
  
  my @files = glob("$indir/s_$lane_nr\_1_*_qseq.txt");

  my ($in1, $out1, $in2, $out2) =(0,0,0,0);

  my ($fh1, $fh2);
  
  $fhs{"$sample_name.1"} = open_outfile( "$basename.1.fq.gz" ) if (! $fhs{"$sample_name.1"} );
  $fhs{"$sample_name.2"} = open_outfile( "$basename.2.fq.gz" ) if (! $fhs{"$sample_name.2"} && -e "$indir/s_$lane_nr\_3_0001_qseq.txt");
  $fhs{"$sample_name.2"} = open_outfile( "$basename.2.fq.gz" ) if (! $fhs{"$sample_name.2"} && ( -e "$indir/s_$lane_nr\_2_0001_qseq.txt" && ! $indexed_run));
  
  foreach my $file (@files) {
    my ($ti, $to) = analyse_tile( $file, $fhs{"$sample_name.1"}, undef, $sample_names{$lane_nr} );
    $in1  += $ti;
    $out1 += $to;
    $file =~ s/(s_\d)_1_/$1_3_/ if ( $indexed_run );
    $file =~ s/(s_\d)_1_/$1_2_/ if ( !$indexed_run );
    if ( -e $file ) {
      # find the next file

      my ($ti, $to) = analyse_tile( $file, $fhs{"$sample_name.2"}, undef, $sample_names{$lane_nr} );
      $in2  += $ti;
      $out2 += $to;
    }

    last if ($debug);
  }


  if($datamonger) {
    EASIH::DONE::add_illumina_lane_stats( $rid, $fids{"$sample_name.1"}, $lane_nr, 1, $sample_name, $in1, $out1 );
    EASIH::DONE::add_illumina_lane_stats( $rid, $fids{"$sample_name.2"}, $lane_nr, 2, $sample_name, $in2, $out2 ) if($in2);

    $reads_pr_sample{$fids{"$sample_name.1"}} += $out1;
    $reads_pr_sample{$fids{"$sample_name.2"}} += $out2 if ( $out2 );
  }
					    
  printf("lane $lane_nr.1\t$sample_name\t$in1\t$out1 (%.2f %%)\t%.2f avg clusters per tile\n", $out1*100/$in1, $out1/120) ;
  printf("lane $lane_nr.2\t$sample_name\t$in1\t$out1 (%.2f %%)\t%.2f avg clusters per tile\n", $out2*100/$in2, $out2/120) if($in2);

  return ($in1, $out1, $in2, $out2);
}


# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub analyse_barcoded_lane {
  my ( $lane_nr) = @_;

  my @files = glob("$indir/s_$lane_nr\_2_*_qseq.txt");

  my %multiplex_stats;
  
  my %barcodes;
  foreach my $lane (sort keys %sample_names) {
    next if (ref ($sample_names{$lane}) ne "HASH");

    foreach my $bcode (keys %{$sample_names{$lane}}) {
#      print "$lane - $bcode ==  $sample_names{$lane}{$bcode} \n";
      my $sample_name = $sample_names{$lane_nr}{$bcode};
      my $basename = $sample_names{ $sample_name };

      $fhs{"$sample_name.1"} = open_outfile( "$basename.1.fq.gz" ) if (! $fhs{"$sample_name.1"} );
      $fhs{"$sample_name.2"} = open_outfile( "$basename.2.fq.gz" ) if (! $fhs{"$sample_name.2"} && -e "$indir/s_$lane_nr\_3_0001_qseq.txt");

    }
  }
  
  my ($in1, $out1, $in2, $out2) =(0,0,0,0);

  foreach my $file (@files) {
    my ($demultiplexing, $counts) = demultiplex_tile($file, \%barcodes);

    map { $multiplex_stats{pass}{$_} += $$counts{pass}{ $_}||0} keys %{$$counts{pass}};
    map { $multiplex_stats{fail}{$_} += $$counts{fail}{ $_}||0} keys %{$$counts{fail}};
    $multiplex_stats{total} += $$counts{total} || 0;

    $file =~ s/(s_\d)_2_/$1_1_/;

    my ($ti, $to) = analyse_tile( $file, undef, $demultiplexing, $sample_names{$lane_nr} );
    $in1  += $ti;
    $out1 += $to;
    $file =~ s/(s_\d)_1_/$1_3_/;
    if ( -e $file  ) {
      # find the next file

      my ($ti, $to) = analyse_tile( $file, undef, $demultiplexing, $sample_names{$lane_nr} );
      $in2  += $ti;
      $out2 += $to;
    }

    last if ($debug);
  }


  printf("lane $lane_nr.1\tMULTIPLEXED\t$in1\t$out1 (%.2f %%)\t%.2f avg clusters per tile\n", $out1*100/$in1, $out1/120) ;
  printf("lane $lane_nr.2\tMULTIPLEXED\t$in1\t$out1 (%.2f %%)\t%.2f avg clusters per tile\n", $out2*100/$in2, $out2/120) if($in2);

  if($datamonger) {
    EASIH::DONE::add_illumina_lane_stats( $rid, undef, $lane_nr, 1, "MULTIPLEXED", $in1, $out1 );
    EASIH::DONE::add_illumina_lane_stats( $rid, undef, $lane_nr, 2, "MULTIPLEXED", $in2, $out2 ) if($in2);
  }

#  print Dumper(\%multiplex_stats);

  foreach my $k ( keys %{$multiplex_stats{pass}} ) {
    my $sample_name = $sample_names{ $lane_nr }{ $k };
    my $perc = sprintf("%.2f", $multiplex_stats{pass}{$k}*100/$multiplex_stats{total});
    printf("lane $lane_nr\t$sample_name\t$k\t$multiplex_stats{pass}{$k}\t$perc %%\n");

    if ($datamonger) {
      EASIH::DONE::add_illumina_multiplex_stats( $rid, $fids{"$sample_name.1"}, $lane_nr, 1, $sample_name, $k, $multiplex_stats{pass}{$k}, $multiplex_stats{pass}{$k} + ($multiplex_stats{fail}{$k} || 0), $perc);
      EASIH::DONE::add_illumina_multiplex_stats( $rid, $fids{"$sample_name.1"}, $lane_nr, 2, $sample_name, $k, $multiplex_stats{pass}{$k}, $multiplex_stats{pass}{$k} + ($multiplex_stats{fail}{$k} || 0), $perc) if ($in2);
    }


    $reads_pr_sample{$fids{"$sample_name.1"}} += $multiplex_stats{pass}{$k} if ( $datamonger);
    $reads_pr_sample{$fids{"$sample_name.2"}} += $multiplex_stats{pass}{$k} if ( $datamonger);

  }

}


# 
# 
# 
# Kim Brugger (09 Jun 2011)
sub analyse_tile {
  my ($input_file, $fout, $demultiplexing, $ebcs) = @_;

  if ( ref($ebcs) eq "HASH") {
    my %found_ebcs;
    foreach my $ebc ( keys %$ebcs ) {
      $found_ebcs{ $ebc } = $$ebcs{$ebc} if ($ebc !~ /^[ACGT]+\z/);
    }
    $ebcs = \%found_ebcs;

    $ebcs = undef if ( keys %$ebcs == 0);

  }
  else {
    $ebcs = undef;
  }
   
  my ( $count_in, $count_out ) = (0,0,0);

  open (my $in, "$input_file") || fail( "Could not open '$input_file': $!\n", "BASECALL2FQ_PATH_ERROR");

  while (my $line = <$in>) {
    chomp $line;

    my @read;
    $count_in++;

    my ($instr, $run_id, $lane, $tile, $x, $y, $index, $read, $bases, $q_line, $filter) = split /\t/, $line;

    $read = 2 if ($read == 3 and $indexed_run);

    #Did not pass the chastity filter
    next if ( $filter == 0);
	
    $bases =~ tr/./N/;           # turn dots into Ns
    $q_line =~ tr/!-\175/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-\136/;

    push @read, "\@${instr}_$run_id:$lane:$tile:$x:$y/$read\n";
    push @read, "$bases\n";
    push @read, "+\n";
    push @read, "$q_line\n";


    if ( $ebcs ) {


      my ( $ebarcode, $decoded, $ebases, $eq_line ) = EASIH::Barcodes::decode_m13f($bases, $q_line);

      # Found the tag, but cannot resolve the barcode
      goto PRINTED if ( $decoded == 0);

      if ( $decoded >= 1 ) {
	print "$ebarcode, $decoded, $ebases, $eq_line\n";

	my $sample_name = $$ebcs { $ebc_set };
	$sample_name .= "_$ebarcode.$read";

	if (! $fhs{"$sample_name"} ) {
	  my $root_dir = "/data/";
	  my $project = substr($sample_name, 0, 3);
	  my $outfile = "$root_dir/$project/raw/$sample_name.fq.gz";
	  if ( $outdir ) {
	    $root_dir = $outdir;
	    $outfile   = "$outdir/$sample_name.fq.gz";
	  }
	  
	  print "$outfile\n";
	  
	  $fhs{"$sample_name"} = open_outfile( "$outfile" );
	}

	my $fh = $fhs{"$sample_name"};
	
	my @read = ("\@${instr}_$run_id:$lane:$tile:$x:$y/$read\n", 
		    "$ebases\n",
		    "+\n",
		    "$eq_line\n");
	
	print $fh join("", @read);
	
	goto PRINTED;
      }
    }

    if ($demultiplexing) {
      my $barcode = $$demultiplexing{ "\@${instr}_$run_id:$lane:$tile:$x:$y"};
      next if ( !$barcode );
      my $sample_name = $sample_names{$lane}{ $barcode };

      my $bfout = $fhs{ "$sample_name.$read" };
      if ( ! $bfout ) {
	print Dumper( \%fhs );
      }
      print $bfout join("", @read);
    }
    else {
      print $fout join("", @read);
    }
  PRINTED:
    $count_out++;
  }
  
  return ($count_in, $count_out);
}



# 
# 
# 
# Kim Brugger (09 Jun 2011)
sub readin_sample_sheet {
  my ( $sample_sheet) = @_;
  
  my ($res, $errors ) = EASIH::Illumina::Sample_sheet::readin( $sample_sheet );

  fail( $errors, "MALFORMED_SAMPLESHEET" ) if ($errors);

  return %$res;
}


# 
# 
# 
# Kim Brugger (10 Aug 2010)
sub validate_lane_names {
  my (%sample_names) = @_;

  if ($sample_sheet ) {
    my $errors = EASIH::Illumina::Sample_sheet::validate( \%sample_names, $limited_lanes ) ;
    fail( $errors, "MALFORMED_SAMPLESHEET") if ( $errors );
  }


  my %basenames;
  for ( my $lane =1; $lane <=8;$lane++) {
    
    next if ( ! $sample_names{ $lane });

    if (ref ($sample_names{$lane}) eq "HASH") {
      foreach my $bcode (keys %{$sample_names{$lane}}) {
	$basenames{ $sample_names{$lane}{$bcode}} = -1;
      }
    }
    else {
      $basenames{ $sample_names{$lane} } = -1;
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

 
  foreach my $basename ( keys %basenames ) {        

    my $root_dir = "/data/";
    my $project = substr($basename, 0, 3);
    my $file = "$root_dir/$project/raw/$basename";

    if ( $outdir ) {
      $root_dir = $outdir;
      $file     = "$outdir/$basename";
    }
    else {
      if ( -e "$root_dir/$project" && ! -d "$root_dir/$project") {
	fail("$root_dir/$project is not a directory\n", "BASECALL2FQ_PATH_ERROR");
      }
      elsif ( -e "$root_dir/$project" && ! -w "$root_dir/$project") {
	fail( "$root_dir/$project is not writeable\n", "BASECALL2FQ_PATH_ERROR");
      }
      
      $root_dir .= "$project/raw/";
      system "mkdir -p $root_dir" if ( ! -d "$root_dir" );
    }
    

    my @files = `find $root_dir | grep $basename | grep fq`;
    my $version = 0;
    if ( @files ) {
      while ( $_ = pop @files ) {
	chomp;

	if ( $version == 0 &&  /[A-Z][0-9]{6,7}.\d+.fq/ ) {
	  $version = 1;
	}
	elsif (  /[A-Z][0-9]{6,7}\_(\d+).\d+.fq/ ) {
	  $version = $1+1 if ($version <= $1 );
	}
	
      }

      $file = "$root_dir/$basename\_$version";
    }
    $sample_names{$basename} = $file;
  }

#  print Dumper( \%sample_names );
#  sleep 60;

  return %sample_names;
}




# 
# As demultiplexing sucks, this have to be done on a tile basis otherwise the memory 
# usage is going to be silly.
# 
# Kim Brugger (04 Jan 2011)
sub demultiplex_tile {
  my ( $file ) = @_;
  
  my ($in1, $out1, $notmplexed1)  = (0, 0, 0);
  my ($in2, $out2, $notmplexed2)  = (0, 0, 0);

  my (%res, %counts);
  my @codes;

  open (my $input, "$file") || fail( "Could not open '$file': $!\n", "BASECALL2FQ_PATH_ERROR");
  while (my $line = <$input>) {

    my ($instr, $run_id, $lane, $tile, $x, $y, $index, $read, $bc, $q_line, $filter) = split /\t/, $line;
    
    #Did not pass the chastity filter
    next if ( $filter == 0);
    # fetch all the barcodes for this lane, once.
    @codes = sort keys %{$sample_names{ $lane }} if (! @codes );
    chop($bc);
    $bc = verify_bcode($bc, @codes);
    
    $counts{total}++ if ($filter);
    
    if ( $bc ) {
      $res{ "\@${instr}_$run_id:$lane:$tile:$x:$y" } = $bc;
      $counts{pass}{$bc}++ if (  $filter );
      $counts{fail}{$bc}++ if ( !$filter );
    }
  }
  close ($input);
  
  return (\%res, \%counts);
}



# 
# 
# 
# Kim Brugger (06 Jan 2011)
sub verify_bcode {
  my ($bc1, @bc2s) = @_;

#  print "$bc1 -- @bc2s\n";

  foreach my $bc2 ( @bc2s ) {

    if ( $no_mismatches) {
      return $bc1 if ($bc1 eq $bc2);
      next;
    }
    else {
      my @seq1 = split('', $bc1);
      my @seq2 = split('', $bc2);
      
      my $diffs = 0;
      for( my $i = 0; $i < @seq1; $i++) {
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

