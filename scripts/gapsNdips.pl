#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (21 Jul 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;


my %opts;
getopts('b:m:hg:l:L:B:', \%opts);

usage() if ( $opts{h});

my $bam_file  = $opts{b} || usage();
my $min_depth = 10;
$min_depth = $opts{m} if ( defined $opts{m});
my $gap_file  = $opts{g};
my $low_file  = $opts{l};
my $region    = $opts{L};
my $regions   = readin_bed( $opts{B} ) if ( $opts{B} );

my ($gout, $lout);

$gap_file = undef if ( $gap_file && $low_file &&  $gap_file eq $low_file);

open ( $gout, "> $gap_file" ) || die "Could not open file '$gap_file': $!\n" if ( $gap_file);
open ( $lout, "> $low_file" ) || die "Could not open file '$low_file': $!\n" if ( $low_file);

my $samtools  = find_program('samtools');
my $bam2depth = find_program('bam2depth');

my $seqs = names_n_lengths( $bam_file );

my (@gaps, @low_coverage);

if ( $regions ) {
  foreach my $name ( keys %$regions ) {
    foreach my $se ( @{$$regions{$name}}) {    
      $name =~ s/chr//i;
      $name ="chr$name";
      my ($start, $end) = @$se;
      my $region = "$name:$start-$end";
      open (my $depth, "$bam2depth $bam_file $region | ") || die "Could not open bam2depth pipeline with $bam_file $name: $!\n";
      analyse_depth( $depth, $name, $start, $end);
    }

    print_report(\@gaps, $gout);
    print_report(\@low_coverage, $lout);
    @gaps = undef;
    @low_coverage = undef;
  }
}
elsif ( $region ) {

  my ($name, $start, $end) = $region =~ /^(\w+?):(\d+)-(\d+)/;
  die "region should be in name:start-end format not '$region'\n"
      if ( ! $name && $start && $end);

  open (my $depth, "$bam2depth $bam_file $region | ") || die "Could not open bam2depth pipeline with $bam_file $name: $!\n";
  analyse_depth( $depth, $name, $start, $end);


}
else {
  foreach my $seq ( @$seqs ) {
  
    my ($name, $length) = @$seq;
    open (my $depth, "$bam2depth $bam_file $name | ") || die "Could not open bam2depth pipeline with $bam_file $name: $!\n";
    analyse_depth( $depth, $name, 0, $length );

  }
}

print_report(\@gaps, $gout);
print_report(\@low_coverage, $lout);



# 
# 
# 
# Kim Brugger (11 Aug 2010)
sub analyse_depth {
  my ($fh, $name, $start, $end) = @_;

#  print "Analysing $name:$start-$end\n";

  my $pre_end = $start;
  my ($low_start, $low_end, $low_depth) = (-1, undef, 0);
  while (<$fh>) {
    chomp;
    my( $region, $start, $depth) = split("\t", $_);
      
    push @gaps, [$name, $pre_end + 1, $start - 1] if ( $start > $pre_end + 1 );
      
    $low_start = $start  if ( $depth < $min_depth && $low_start == -1);
    $low_end = $start    if ( $depth < $min_depth);
    $low_depth += $depth if ( $depth < $min_depth);
    
    if ( $low_start != -1 && $depth > $min_depth ) {
      push @low_coverage, [$name, $low_start, $start -1, int($low_depth/($start - $low_start))];
      $low_start = -1;
      $low_depth = 0;
    }
    
    $pre_end = $start;
  }  
  
  push @low_coverage, [$name, $low_start, $low_end, int($low_depth/($low_end - $low_start + 1))] if ( $low_start != -1 );
  
  push @gaps, [$name, $pre_end + 1, $end] if ($pre_end + 1 < $end );
}




# 
# 
# 
# Kim Brugger (21 Jul 2010)
sub print_report {
  my ($entries, $stream ) = @_;

  $stream = *STDOUT if ( ! $stream );

  foreach my $entry ( @$entries ) {
    next if (!$entry);
    $$entry[3] ||= 0;
    print $stream "$$entry[0]:$$entry[1]-$$entry[2]\t$$entry[3]\n";
#    print $stream join("\t",@$entry)."\n";
  }
}





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

    chomp;
    my ($chr, $start, $end) = split("\t", $_);

    ($chr, $start, $end) = $_ =~ /(.*?):(\d+)-(\d+)/
	if ( ! $start );

    next if ( ! $chr);
    
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
