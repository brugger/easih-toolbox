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

my %opts;
getopts("a:1:2:3:4:5:6:7:8:hs:i:o:Sc:lhb:M", \%opts);

usage() if ( ! $opts{a} && ! $opts{1} && ! $opts{2} && ! $opts{3} && ! $opts{4} && 
	   ! $opts{5} && ! $opts{6} && ! $opts{7} && ! $opts{8}  ||  $opts{h});

my $parallel =  1;
$parallel = 0 if ($opts{S});
my $limited_lanes = $opts{'l'};
my $MAX_NODES = $opts{c} || 8;
my $allow_mismatches = $opts{M};
my %barcoded_lanes;
map { $barcoded_lanes{$_}++} split('', $opts{b}) if ($opts{b});

# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub usage {
  $0 =~ s/.*\///;
  print "USAGE: $0 extracts data from the illumina Base directory\n";
  print "USAGE: -1,  -2, .., -8 < the lanes are assigned this name>\n";
  print "USAGE: -a< id, all lanes that are not named specifically gets this id>\n";  
  print "USAGE: -S<erial version>\n";
  print "USAGE: -c[pus to run on if runing in parallel, default 7]\n";
  print "USAGE: -l<imited lanes, by default the whole slide is extracted>\n";
  print "USAGE: -h<elp>\n";
  print "USAGE: -b[arcoded lanes eg: -b 345 for lane 3, 4, 5]\n";
  exit -1;
}


my $lane1 = $opts{'1'};
my $lane2 = $opts{'2'};
my $lane3 = $opts{'3'};
my $lane4 = $opts{'4'};
my $lane5 = $opts{'5'};
my $lane6 = $opts{'6'};
my $lane7 = $opts{'7'};
my $lane8 = $opts{'8'};

$lane1 = $opts{a} if ($opts{a} && !$opts{'1'});
$lane2 = $opts{a} if ($opts{a} && !$opts{'2'});
$lane3 = $opts{a} if ($opts{a} && !$opts{'3'});
$lane4 = $opts{a} if ($opts{a} && !$opts{'4'});
$lane5 = $opts{a} if ($opts{a} && !$opts{'5'});
$lane6 = $opts{a} if ($opts{a} && !$opts{'6'});
$lane7 = $opts{a} if ($opts{a} && !$opts{'7'});
$lane8 = $opts{a} if ($opts{a} && !$opts{'8'});

my @lanes = (undef, $lane1, $lane2, $lane3, $lane4, $lane5, $lane6, $lane7, $lane8);
my $indir   = $opts{'i'} || "./";
my $outdir  = $opts{'o'} || "./";
system "mkdir $outdir" if ( ! -d $outdir);

my @fhs       = ();
my $total     = 0;

validate_lane_names();

my @params;

for(my $i = 1; $i<=8; $i++) {

  my $lane_name = $lanes[ $i ];
  
  next if ( ! $lane_name );

  if ( $opts{b} ) {
    # as the illumina is stupid the barcode is in lane two.

    my $dm;
    if ($barcoded_lanes{$i}) {
      push @params, [\&demultiplex, $i, $lane_name];
    }
    else {
      push @params, [\&analyse_lane, $i, "1", $lane_name, "$indir/s_$i\_1_*_qseq.txt"];
      push @params, [\&analyse_lane, $i, "2", $lane_name, "$indir/s_$i\_3_*_qseq.txt"];
    }
  }
  # This is standard, how things normally are...
  else {
    push @params, [\&analyse_lane, $i, "1", $lane_name, "$indir/s_$i\_1_*_qseq.txt"];
    push @params, [\&analyse_lane, $i, "2", $lane_name, "$indir/s_$i\_2_*_qseq.txt"];
  }
}

if ( $parallel ) {
  run_parallel( @params);
}
else {
  run_serial( @params );
}


# 
# As demultiplexing sucks, this have to be done on a tile basis otherwise the memory 
# usage is going to be silly.
# 
# Kim Brugger (04 Jan 2011)
sub demultiplex {
  my ($lane, $lane_name) = @_;
  
  my @files = glob("$indir/s_$lane\_2_*_qseq.txt");

  my ($in1, $out1, $notmplexed1)  = (0, 0, 0);
  my ($in2, $out2, $notmplexed2)  = (0, 0, 0);

  foreach my $file (  @files) {
    open (my $input, "$file") || die "Could not open '$file': $!\n";
    my %res;
    
    while (my $line = <$input>) {

      my ($instr, $run_id, $lane, $tile, $x, $y, $index, $read, $bc, $q_line, $filter)
	  = split /\t/, $line;

      #Did not pass the chastity filter
      next if ( $filter == 0);
      push @{$res{$bc}}, "\@${instr}_$run_id:$lane:$tile:$x:$y";
    }
    close ($input);

    my (%good_bc, %codes);
    foreach my $bc ( sort {@{$res{$b}} <=> @{$res{$a}}} keys %res) {

      next if ( $bc =~ /^\.+\Z/);
      
      if ( ! keys %good_bc ) {
	map {$good_bc{ $_ } = $bc } @{$res{ $bc }};
	$codes{ $bc }++;
      }
      else {
	# note to self: this is the ugliest code I have produced in ages. Well done, go pamper yourself!
	if ( $allow_mismatches ) {
	  my $bc2 = similar_bcs( $bc, keys %codes);
	  if ( $bc2 ) {
	    map {$good_bc{ $_ } = $bc2 } @{$res{ $bc }};
	    $codes{ $bc2 }++;
	    next;
	  }
	}

	next if ( $bc =~ /\./);
	 
	if ( @{$res{$bc}} > 500 ) {
	  map {$good_bc{ $_ } = $bc } @{$res{ $bc }};
	  $codes{ $bc }++;
	}
      }
    }


    $file =~ s/_2_/_1_/;
    my ($in, $out, $notmplexed) = analyse_lane($lane, "1", $lane_name, $file, \%good_bc, [keys %codes] );
    $in1  += $in;
    $out1 += $out;
    $notmplexed1 += $notmplexed;
    $file =~ s/_1_/_3_/;
    ($in, $out, $notmplexed) = analyse_lane($lane, "2", $lane_name, $file, \%good_bc, [keys %codes] );
    $in2  += $in;
    $out2 += $out;
    $notmplexed2 += $notmplexed;
  }

  printf STDERR ("lane $lane.1\t$lane_name\t$in1\t$out1 (%.2f %%)\t%.2f avg clusters per tile", $out1*100/$in1, $out1/120) ;
  printf STDERR ("\tnot demultiplexed: %d (%.2f %%)\n", $notmplexed1, $notmplexed1*100/$in1);

  if ( $out2) {
    printf STDOUT ("lane $lane.2\t$lane_name\t$in2\t$out2 (%.2f %%)\t%.2f avg clusters per tile", $out2*100/$in2, $out2/120) ;
    printf STDOUT ("\tnot demultiplexed: %d (%.2f %%)\n", $notmplexed2, $notmplexed2*100/$in2);
  }

}


# 
# 
# 
# Kim Brugger (06 Jan 2011)
sub similar_bcs {
  my ($bc1, @bc2s) = @_;

  foreach my $bc2 ( @bc2s ) {
    my @seq1 = split('', $bc1);
    my @seq2 = split('', $bc2);

    my $diffs = 0;
    for( my $i = 0; $i < @seq1; $i++) {
      $diffs++ if ( $seq1[$i] ne $seq2[$i]);
    }

    return $bc2 if ( $diffs <= 1);
  }

  return undef;
}



# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub run_serial {
  my (@params) = @_;

  while( my $param = shift @params) {
    my $function = shift @$param;
    &$function( @$param);
  }
}



# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub run_parallel {
  my (@params) = @_;

  use POSIX ':sys_wait_h';
  my @cpids     = ();
  
  my $done = 0; # to track the number of files handled

  my $running_nodes = 0;

  while (1) {
  FREE_NODE:
    if ($running_nodes <  $MAX_NODES) {
      my $param = shift @params;

      last if ( !$param);

      my $function = shift @$param;
      
      my $cpid = create_child($function, @$param);
      $running_nodes++;
      $total++;
      push @cpids, $cpid;
    }
    else {
      # loop through the nodes to see when one becomes available ...
      while ($done < $total) {
	for (my $i = 0; $i <@cpids; $i++) {
	  next if ($cpids[$i] == -10);
	  
	  my $cpid = $cpids[$i];
	  if (!waitpid($$cpid, WNOHANG)) {
#	  print "Waiting for ($$cpid)\n";
	  }
	  elsif ($$cpid != -10) {
	    $done++;
	    $cpids[$i] = -10;
	    $running_nodes--;
	  }
	}
	sleep 10;
	
	last if ($running_nodes < $MAX_NODES);
      }
      goto FREE_NODE;
    }
    
    last if ($running_nodes == 0);
  }

  while ($done < $total) {
    for (my $i = 0; $i <@cpids; $i++) {
      next if ($cpids[$i] == -10);
    
      my $cpid = $cpids[$i];
      if (!waitpid($$cpid, WNOHANG)) {
	;
      }
      elsif ($$cpid != -10) {
	$done++;
	$cpids[$i] = -10;
      }
    }
    sleep 10;
  }

  foreach my $fh (@fhs ) {
    print <$fh>;
    close $fh;
  }

}



sub create_child {
  my ($command, @params) = @_;

  my $pid;
  # open a pipe from the kid to the parent so they can communicate
  if ($pid = open($fhs[$total], "-|")) {
    ;
  } 
  else {
    die "cannot fork: $!" unless defined $pid;
    &$command(@params);
    exit;
  }
  
  return \$pid;
}




# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub analyse_lane {
  my ( $lane_nr, $lane_index, $lane_name, $glob, $multiplex, $barcodes) = @_;

  my @files = glob($glob);
  my ($count_in, $count_out, $not_demultiplexed) = (0,0,0);
  my (%fhs, $out);
  if ( @files ) {
    if ( $barcodes ) {
      foreach my $barcode ( @$barcodes) {
	open ($fhs{$barcode}, "| gzip -c >> $outdir/$lane_name\_$barcode.$lane_index.fq.gz") || die "Could not open '$outdir/$lane_name\_$barcode.$lane_index.fq.gz': $!\n";
      }
    }
    else {
      open ($out, "| gzip -c > $outdir/$lane_name.$lane_index.fq.gz") || die "Could not open '$outdir/$lane_name.$lane_index.fq.gz': $!\n";
    }


    foreach my $file (  @files) {
      open (my $in, "$file") || die "Could not open '$file': $!\n";

      my (@reads, %bc_reads);
    
      while (my $line = <$in>) {

	my @read;
	chomp $line;
	$count_in++;

	my ($instr, $run_id, $lane, $tile, $x, $y, $index, $read, $bases, $q_line, $filter) = split /\t/, $line;

	#Did not pass the chastity filter
	next if ( $filter == 0);
	
	$bases =~ tr/./N/;           # turn dots into Ns
	$q_line =~ tr/!-\175/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-\136/;

	if ($index ne '0') {
	  push @read, "\@${instr}_$run_id:$lane:$tile:$x:$y\#$index/$read\n";
	} else {
	  push @read, "\@${instr}_$run_id:$lane:$tile:$x:$y/$read\n";
	}
	push @read, "$bases\n";
	push @read, "+\n";
	push @read, "$q_line\n";
	
	if ($multiplex) {
	  my $barcode = $$multiplex{ "\@${instr}_$run_id:$lane:$tile:$x:$y"};
	  if ( !$barcode ) {
	    $not_demultiplexed++;
	    next;
	  }
	  $count_out++;
	  push @{$bc_reads{$barcode}}, \@read;
	}
	else {
	  $count_out++;
	  push @reads, \@read;
	}
      }
      
      if ( $multiplex ) {
	foreach my $barcode ( keys %bc_reads ) {
	  foreach my $read (sort {$$a[0] cmp $$b[0]} @{$bc_reads{ $barcode }}) {
	    print {$fhs{$barcode}} join("", @$read);
	  }
	}
      }
      else {
	foreach my $read (sort {$$a[0] cmp $$b[0]} @reads) {
	  print $out join("", @$read);
	}
      }
    }
    
    if ( ! $multiplex ) {
      printf STDOUT ("lane $lane_nr.$lane_index\t$lane_name\t$count_in\t$count_out (%.2f %%)\t%.2f avg clusters per tile", $count_out*100/$count_in, $count_out/120) ;
      print  STDOUT ("\n");
    }
    else {
      printf STDOUT ("lane $lane_nr.$lane_index\t$lane_name\t$count_in\t$count_out (%.2f %%)\t%.2f avg clusters per tile", $count_out*100/$count_in, $count_out/120) ;
      printf STDOUT ("\tnot demultiplexed: %d (%.2f %%)", $not_demultiplexed, $not_demultiplexed*100/$count_in) if ($multiplex);
      print  STDOUT ("\n");
    }


  }
  close ($out) if ( $out);

  if ( $barcodes ) {
    foreach my $barcode ( @$barcodes) {
      close ($fhs{$barcode});
    }
  }

  
  return ($count_in, $count_out, $not_demultiplexed);
}



# 
# 
# 
# Kim Brugger (10 Aug 2010)
sub validate_lane_names {

  my (%seen, $error);
  my @missing_names;
  for( my $i=1; $i< @lanes; $i++) {
    my $name = $lanes[ $i ];
    if ( ! defined $name || $name eq "") {
      push @missing_names, $i;
      next;
    }
    push @{$seen{$name}}, $i;
  }

  die "Missing names for the following lanes: @missing_names\n" if ( @missing_names && ! $limited_lanes);

 
  foreach my $name ( keys %seen) {
    
    if ( @{$seen{$name}} > 1)  {
      my @new_names;
      my $i = 1;
      map {$lanes[ $_ ] .= "_$i"; push @new_names, $lanes[ $i++ ]} @{$seen{$name}};
      print "$name is being used multiple times. Renaming $name to @new_names \n";
    }
  }
}

