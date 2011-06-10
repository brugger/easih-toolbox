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

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::Git;

my %opts;
getopts("a:1:2:3:4:5:6:7:8:hs:i:o:lhn", \%opts);


my $limited_lanes = $opts{'l'};
my $no_mismatches = $opts{n};

# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub usage {
  $0 =~ s/.*\///;
  print "USAGE: $0 extracts data from the illumina Base directory\n";
  print "USAGE: -1,  -2, .., -8 < the lanes are assigned this name> Overrules the sample_sheet.csv\n";
  print "USAGE: -a< id, all lanes that are not named specifically gets this id>\n";  
  print "USAGE: -h<elp>\n";
  print "USAGE: -i<nput dir, default BaseCalls>\n";
  print "USAGE: -l<imited lanes, by default the whole slide is extracted>\n";
  print "USAGE: -n<o mismatches in barcodes, normal is 1 error>\n";
  print "USAGE: -o<utput dir, default is /data/<project ID>/raw/, on a sample basis>\n";
  print "USAGE: for barcoded lanes, please use a sample sheet\n";
  exit -1;
}

my $indir       = $opts{'i'} || "./";
my $outdir      = $opts{'o'};

my $indexed_run = 0;
my $sample_sheet = $opts{'s'};
$sample_sheet = "$indir/sample_sheet.csv" if (!$sample_sheet && -e "$indir/sample_sheet.csv");
if (!$sample_sheet && -e "BaseCalls/sample_sheet.csv") {
  $indir = "BaseCalls";
  $sample_sheet = "$indir/sample_sheet.csv";
}
if (!$sample_sheet && -e "sample_sheet.csv") {
  $indir = "./";
  $sample_sheet = "sample_sheet.csv";
}

usage() if (! $sample_sheet &&  ! $opts{a} && ! $opts{1} && ! $opts{2} && ! $opts{3} && ! $opts{4} && 
	   ! $opts{5} && ! $opts{6} && ! $opts{7} && ! $opts{8}  ||  $opts{h});



my %sample_names = readin_sample_sheet( $sample_sheet);

$sample_names{1} = $opts{'1'} if ($opts{'1'});
$sample_names{2} = $opts{'2'} if ($opts{'2'});
$sample_names{3} = $opts{'3'} if ($opts{'3'});
$sample_names{4} = $opts{'4'} if ($opts{'4'});
$sample_names{5} = $opts{'5'} if ($opts{'5'});
$sample_names{6} = $opts{'6'} if ($opts{'6'});
$sample_names{7} = $opts{'7'} if ($opts{'7'});
$sample_names{8} = $opts{'8'} if ($opts{'8'});

$sample_names{1} = $opts{a} if ($opts{a} && !$opts{'1'});
$sample_names{2} = $opts{a} if ($opts{a} && !$opts{'2'});
$sample_names{3} = $opts{a} if ($opts{a} && !$opts{'3'});
$sample_names{4} = $opts{a} if ($opts{a} && !$opts{'4'});
$sample_names{5} = $opts{a} if ($opts{a} && !$opts{'5'});
$sample_names{6} = $opts{a} if ($opts{a} && !$opts{'6'});
$sample_names{7} = $opts{a} if ($opts{a} && !$opts{'7'});
$sample_names{8} = $opts{a} if ($opts{a} && !$opts{'8'});



%sample_names = validate_lane_names(%sample_names);
my %fhs;

for(my $lane = 1; $lane<=8; $lane++) {
 
  # this lane is barcoded...
  if (ref ($sample_names{$lane}) eq "HASH") {
    analyse_barcoded_lane($lane);
  }
  else {
    analyse_lane($lane)
  }

}

# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub analyse_lane {
  my ( $lane_nr) = @_;

  my $sample_name = $sample_names{$lane_nr};
  my $basename = $sample_names{ $sample_name };
  
  my @files = glob("$indir/s_$lane_nr\_1_*_qseq.txt");

  my ($in1, $out1, $in2, $out2) =(0,0,0,0);

  my ($fh1, $fh2);
  
  open ($fhs{"$sample_name.1"}, "| gzip -c > $basename.1.fq.gz") || die "Could not open '$basename.1.fq.gz': $!\n"
      if (! $fhs{"$sample_name.1"} );

  open ($fhs{"$sample_name.2"}, "| gzip -c > $basename.2.fq.gz") || die "Could not open '$basename.2.fq.gz': $!\n"
      if (! $fhs{"$sample_name.2"} && -e "$indir/s_$lane_nr\_3_0001_qseq.txt");

  open ($fhs{"$sample_name.2"}, "| gzip -c > $basename.2.fq.gz") || die "Could not open '$basename.2.fq.gz': $!\n"
      if (! $fhs{"$sample_name.2"} && ( -e "$indir/s_$lane_nr\_2_0001_qseq.txt" && ! $indexed_run));
  
  foreach my $file (@files) {
    my ($ti, $to) = analyse_tile( $file, $fhs{"$sample_name.1"} );
    $in1  += $ti;
    $out1 += $to;
    $file =~ s/(s_\d)_1_/$1_3_/ if ( $indexed_run );
    $file =~ s/(s_\d)_1_/$1_2_/ if ( !$indexed_run );
    if ( -e $file ) {
      # find the next file

      my ($ti, $to) = analyse_tile( $file, $fhs{"$sample_name.2"} );
      $in2  += $ti;
      $out2 += $to;
    }

    last;
  }

#  printf("lane $lane_nr.1\t$sample_name\t$in1\t$out1\n") ;
#  printf("lane $lane_nr.2\t$sample_name\t$in2\t$out2\n");



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
      open ($fhs{"$sample_name.1"}, "| gzip -c > $basename.1.fq.gz") || die "Could not open '$basename.1.fq.gz': $!\n"
	  if (! $fhs{"$sample_name.1"} );

      open ($fhs{"$sample_name.2"}, "| gzip -c > $basename.2.fq.gz") || die "Could not open '$basename.2.fq.gz': $!\n"
	  if (! $fhs{"$sample_name.2"} && -e "$indir/s_$lane_nr\_3_0001_qseq.txt");

    }
  }
  
  my ($in1, $out1, $in2, $out2) =(0,0,0,0);

  foreach my $file (@files) {
    my ($demultiplexing, $counts) = demultiplex_tile($file, \%barcodes);

    map { $multiplex_stats{$_} += $$counts{ $_}} keys %$counts;

    $file =~ s/(s_\d)_2_/$1_1_/;

    my ($ti, $to) = analyse_tile( $file, undef, $demultiplexing );
    $in1  += $ti;
    $out1 += $to;
    $file =~ s/(s_\d)_1_/$1_3_/;
    if ( -e $file  ) {
      # find the next file

      my ($ti, $to) = analyse_tile( $file, undef, $demultiplexing );
      $in2  += $ti;
      $out2 += $to;
    }

#    last;
  }

  printf("lane $lane_nr.1\tMULTIPLEXED\t$in1\t$out1 (%.2f %%)\t%.2f avg clusters per tile\n", $out1*100/$in1, $out1/120) ;
  printf("lane $lane_nr.2\tMULTIPLEXED\t$in1\t$out1 (%.2f %%)\t%.2f avg clusters per tile\n", $out2*100/$in2, $out2/120) if($in2);

  foreach my $k ( keys %multiplex_stats ) {
    next if ($k eq "total");
    my $sample_name = $sample_names{ $lane_nr }{ $k };
    printf("lane $lane_nr\t$sample_name\t$k\t$multiplex_stats{$k}\t%.2f %%\n", $multiplex_stats{$k}*100/$multiplex_stats{total});
  }

}


# 
# 
# 
# Kim Brugger (09 Jun 2011)
sub analyse_tile {
  my ($input_file, $fout, $demultiplexing) = @_;

  my ( $count_in, $count_out ) = (0,0,0);


  open (my $in, "$input_file") || die "Could not open '$input_file': $!\n";

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
  
  my (%res );
  
  my $text_delim = "";
  my $field_delim = "";

  open(my $in, $sample_sheet) || die "Could not open '$sample_sheet': $!\n";
  my @lines;
  while(<$in>) {
    $_ =~ s/\r\n/\n/g; 
    $_ =~ s/\n\r/\n/g; 
    $_ =~ s/\r/\n/g; 
    push @lines, split("\n",$_);
  }
  close $in;

  while($_ = shift @lines ) {
    chomp;
    
    
    # As I dont trust they can export the csv file in the same format each time
    # we will use the first line to identify field and text delimiters.
    if (/^(.{0,1})FCID/) {
      $text_delim = $1;
      /FCID$text_delim(.)/;
      $field_delim = $1;
    }      
    else {
      my @F = split($field_delim, $_);
      my (undef, $lane, $sample_id, undef, $index, undef) = @F;

      $index ||= "";

      $lane      =~ s/^$text_delim(.*)$text_delim\z/$1/;
      $sample_id =~ s/^$text_delim(.*)$text_delim\z/$1/;
      $index     =~ s/^$text_delim(.*)$text_delim\z/$1/;

      if ( $index ) {
	$res{$lane}{$index} = $sample_id;
	$indexed_run++;
      }
      else {
	$res{$lane} = $sample_id;
      }
    }
  }


  return %res;
}


# 
# 
# 
# Kim Brugger (10 Aug 2010)
sub validate_lane_names {
  my (%sample_names) = @_;

  my %basenames;
  foreach my $lane ( keys %sample_names) {
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
      die "$outdir is not a directory\n";
    }
    elsif ( -e "$outdir" && ! -w "$outdir") {
      die "$outdir is not writeable\n";
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
	die "$root_dir/$project is not a directory\n";
      }
      elsif ( -e "$root_dir/$project" && ! -w "$root_dir/$project") {
	die "$root_dir/$project is not writeable\n";
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

  open (my $input, "$file") || die "Could not open '$file': $!\n";
  while (my $line = <$input>) {

    my ($instr, $run_id, $lane, $tile, $x, $y, $index, $read, $bc, $q_line, $filter) = split /\t/, $line;
    
    #Did not pass the chastity filter
    next if ( $filter == 0);
    # fetch all the barcodes for this lane, once.
    @codes = sort keys %{$sample_names{ $lane }} if (! @codes );
    chop($bc);
    $bc = verify_bcode($bc, @codes);
    
    $counts{total}++;
    
    if ( $bc ) {
      $res{ "\@${instr}_$run_id:$lane:$tile:$x:$y" } = $bc;
      $counts{$bc}++;
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
