#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Nov 2010), contact: kim.brugger@easih.ac.uk
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts('hB:b:', \%opts);

usage() if ( $opts{h});

my $regions   =  $opts{B} || usage();
my $bamfile   =  $opts{b} || usage();

$regions   = readin_bed( $regions );

my $samtools  = find_program('samtools');
my $bam2depth = find_program('bam2depth');

my $seqs = names_n_lengths( $bamfile );

my ($START, $END) = (0, 1);

my $total_reads;
my %base_coverage;
my %exon_coverage;

my $escape = 1;

my $patched_start = 0;

foreach my $chr ( sort keys %$regions ) {

#  $chr ="chr2";

  my (undef, $base_pos, $depth);
  my @regions =  @{$$regions{$chr}};
  my $region  = shift @regions;
  my ($start, $end) = @{$region};

#  print "Doing $chr ... expecting ".@regions." regions\n";

  open (my $pipe, "$bam2depth $bamfile $chr | " ) || die "Could not open pipe: $!\n";
  my $covered = 0;
  while(<$pipe>) {
    chomp;
    (undef, $base_pos, $depth) = split("\t", $_);
    
    next if ( $base_pos < $start);

    if ( $base_pos >= $start && 
	 $base_pos <= $end ) {

      $base_coverage{ $depth }++;
      $covered++;
    }
    else {
      
      if ( $start != $end ) {
      
	$exon_coverage{ int( $covered/($end-$start + 1)*100)}++;
	$base_coverage{ 0 } += $end - $start + 1 - $covered;
      }

      ($start, $end) = (undef, undef);
      $region  = shift @regions;
      last if ( ! $region );
      ($start, $end) = @{$region};

      $covered = 0;
    }      
  }
  close( $pipe );

  if (  $start ) {
      if ( $start != $end ) {
	$exon_coverage{ int( $covered/($end-$start + 1)*100)}++;
	$base_coverage{ 0 } += $end - $start + 1 - $covered;
      }
  }
  
  if (@regions) {
    foreach $region (@regions) {
      ($start, $end) = @{$region};
	$exon_coverage{ 0 }++;
	$base_coverage{ 0 } += $end - $start + 1;
    }
  }

#  last;
}

my %binned = ( 0 => 0,
	       1 => 0,
	       10 => 0,
	       20 => 0,
	       30 => 0,
	       40 => 0);
my $base_total = 0;
foreach my $depth ( keys %base_coverage ) {  

  $base_total += $base_coverage{ $depth };
  if ( $depth == 0 || $depth == 1 ) {
    $binned{ $depth } = $base_coverage{ $depth };
    next;
  }
  elsif ( $depth < 10 ) {
    $binned{ 10 } +=  $base_coverage{ $depth };
  }
  elsif ( $depth < 20 ) {
    $binned{ 20 } +=  $base_coverage{ $depth };
  }
  elsif ( $depth < 30 ) {
    $binned{ 30 } +=  $base_coverage{ $depth };
  }
  else {
    $binned{ 40 } +=  $base_coverage{ $depth };
  }
}

print "$bamfile\n";
print "base depth distribution:\n"; 
printf("     0x depth: $binned{  0 } (%.2f %%) \n", $binned{  0 }/$base_total*100);
printf("     1x depth: $binned{  1 } (%.2f %%) \n", $binned{  1 }/$base_total*100);
printf("  2-10x depth: $binned{ 10 } (%.2f %%) \n", $binned{ 10 }/$base_total*100);
printf(" 10-20x depth: $binned{ 20 } (%.2f %%) \n", $binned{ 20 }/$base_total*100);
printf(" 20-30x depth: $binned{ 30 } (%.2f %%) \n", $binned{ 30 }/$base_total*100);
printf("   >30x depth: $binned{ 40 } (%.2f %%) \n", $binned{ 40 }/$base_total*100);


my %exons = (  0 => 0,
	      70 => 0,
	      80 => 0,
	      90 => 0,
	     100 => 0,
	     101 => 0);

my $exon_total = 0;
foreach my $coverage ( keys %exon_coverage ) {  

  $exon_total += $exon_coverage{ $coverage };
  if ( $coverage == 0 ) {
    $exons{ $coverage } += $exon_coverage{ $coverage };
  }
  elsif ( $coverage < 70 ) {
    $exons{ 70 } += $exon_coverage{ $coverage };
  }
  elsif ( $coverage < 80 ) {
    $exons{ 80 } += $exon_coverage{ $coverage };
  }
  elsif ( $coverage < 90 ) {
    $exons{ 90 } += $exon_coverage{ $coverage };
  }
  elsif ( $coverage < 100 ) {
    $exons{ 100 } += $exon_coverage{ $coverage };
  }
  elsif ( $coverage >= 100 ) {
    $exons{ 101 } += $exon_coverage{ $coverage };
  }
}


print "bait regions coverage:\n"; 
printf("      0 %% of bait covered: $exons{  0 } (%.2f %%) \n", $exons{     0 }/$exon_total*100);
printf("    <70 %% of bait covered: $exons{  70 } (%.2f %%) \n", $exons{   70 }/$exon_total*100);
printf("  70-80 %% of bait covered: $exons{  80 } (%.2f %%) \n", $exons{   80 }/$exon_total*100);
printf("  80-90 %% of bait covered: $exons{  90 } (%.2f %%) \n", $exons{   90 }/$exon_total*100);
printf("  90-99 %% of bait covered: $exons{  100 } (%.2f %%) \n", $exons{  100 }/$exon_total*100);
printf("    100 %% of bait covered: $exons{  101 } (%.2f %%) \n", $exons{  101 }/$exon_total*100);

# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub names_n_lengths {
  my ( $bam_file ) = @_;

  my @sequences = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);

    my ($name, $length);
    foreach my $field ( split("\t") ) {
      $name   = $1 if ( $field =~ /SN:(.*)/);
      $length = $1 if ( $field =~ /LN:(\d+)/);
    }

    push @sequences, [$name, $length] if ( $name && $length );
  }

    
  return \@sequences;
}



# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub usage {
  
  $0 =~ s/.*\///;
  die "USAGE: $0 -b[am file] -m[in depth, default 200] -g[ap out file, default stdout] -l[ow regions file, default stdout] -L[ only look at this region (name:start-end) -B[ed file with regions of interest]\n";
}



# 
# 
# 
# Kim Brugger (13 Jul 2010)
sub find_program {
  my ($program) = @_;


  my @paths = ("/home/easih/bin/",
	       "/home/kb468/bin/",
	       "/home/kb468/easih-toolbox/scripts/",
	       "/usr/local/bin");
  
  foreach my $path ( @paths ) {
    
    return "$path/$program" if ( -e "$path/$program" );
  }

  my $location = `which $program`;
  chomp( $location);
  
  return $location if ( $location );

  return undef;
}


# 
# 
# 
# Kim Brugger (11 May 2010)
sub readin_bed {
  my ( $infile, $merge ) = @_;

  $merge = 1;
  my %res;

  open (STDIN, $infile) || die "Could not open '$infile': $!\n" if ( $infile );
  while(<STDIN>) {

    next if (/^chromosome/);

    chomp;
    my ($chr, $start, $end) = split("\t", $_);

    ($chr, $start, $end) = $_ =~ /(.*?):(\d+)-(\d+)/
	if ( ! $start );

    next if ( ! $chr );

    $chr =~ s/chr//;
#    $chr = "chr$chr";
    
    push @{$res{$chr}}, [$start, $end];
  }

  if ( $merge ) {
    
    foreach my $key ( keys %res ) {
      
      @{$res{$key}} = sort { $$a[0] <=> $$b[0] } @{$res{$key}};
      my @tmp;
      my @data = @{$res{$key}};
      
      for(my $i=0;$i< @data; $i++) {
	
	# need at least one element in the array, so push and move on.
	if ( ! @tmp ) {
	  push @tmp, $data[ $i ];
	  next;
	}
	
	# contained in the region
	if ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ] &&
	     $data[ $i ][ 1 ] <= $tmp[ -1 ][ 1 ]) {
	  next;
	}
	# overlapping
	elsif ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ]  &&
		$data[ $i ][ 0 ] <= $tmp[ -1 ][ 1 ]) {
	  
	  $tmp[ -1 ][ 1 ] = $data[ $i ][ 1 ];
	}
	# There is a gap between the end block and this one. Just push it on the end of  the array!
	else {
	  push @tmp, $data[ $i ];
	}
      }
      @{$res{$key}} = @tmp;
    }
  }

  return \%res;
}
