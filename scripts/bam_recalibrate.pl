#!/usr/bin/perl -w
# 
# Recalibrates QV values by calculating the true error rate for all QVs.
# 
#
# Kim Brugger (11 Oct 2011), contact: kim.brugger@easih.ac.uk

use strict;
use Data::Dumper;
use warnings;
use Getopt::Std;

my %opts;
getopts('b:d:R:B:hs:Sr:UM:m:L:l:s:m:o:', \%opts);

my %calib_stats;
my %calib_matrix;


my $bam_file = $opts{'b'} || Usage();
my $chr_file = $opts{'R'} || Usage();
my $base_mutation_rate = $opts{'m'};
my $filter_on_base_mutation_rate = 0 if ( $base_mutation_rate && 
					  $base_mutation_rate > 0 && $base_mutation_rate < 100);

die "'$chr_file' does not exist\n" 
    if ( ! -e $chr_file );
die "index does not exist for '$chr_file', please create one with samtools faidx\n"
    if (! -e "$chr_file.fai");

my $samtools    = find_program('samtools');

my $MAX_READS = 0;
my $MAX_BASES = 0;

my $outfile        = $opts{'o'};

my $log_file       = $opts{'L'} || 0;
my $report_file    = $opts{'l'};


open (my $log, "> $log_file") || die "Could not open file '$log_file': $!\n" if ( $log_file );

my ( $s, $d, $e, $b_pre, $b_post) = (0,0,0, 0, 0);

my $region      = $opts{'r'};
my $sample_size = $opts{'s'} || 0;
$region = "'gi|170079663|ref|NC_010473.1|'";

my $fasta;

my %hist;

my %QV_stats;

if ( ! $region ) {
  # loop through the regions one by one, but only keep the most current chr in memory.
  
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      if ( $field =~ /SN:(.*)/) {
	$fasta = readfasta( $chr_file, $1 );
	analyse($1)
      }
    }
  }
}
else {
  $region =~ s/,//g;
  if ($region =~ /^(\w+):\d+-\d+/ || $region =~ /^(\w+):\d+\z/ ) {
    $fasta = readfasta( $chr_file, $1 );
    analyse($region)
  }
  else {
    $fasta = readfasta( $chr_file, $region );
    analyse($region)
  }
}


foreach my $QV ( sort {$a <=> $b} keys %calib_stats ) {
#x  print Dumper( $calib_stats{ $QV });
  $calib_matrix{ chr($QV + 33) } = chr(phred_value( $calib_stats{ $QV }{ M }, $calib_stats{ $QV }{ X }) + 33 );
  printf("$QV\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d (%d, %d)\n",
	 phred_value( $calib_stats{ $QV }{ M }, $calib_stats{ $QV }{ X }),
	 phred_value( $calib_stats{ $QV }{ A }{ M }, $calib_stats{ $QV }{ A }{ X }),
	 phred_value( $calib_stats{ $QV }{ C }{ M }, $calib_stats{ $QV }{ C }{ X }),
	 phred_value( $calib_stats{ $QV }{ G }{ M }, $calib_stats{ $QV }{ G }{ X }),
	 phred_value( $calib_stats{ $QV }{ T }{ M }, $calib_stats{ $QV }{ T }{ X }),
	 ($calib_stats{ $QV }{ M } || 0) + ($calib_stats{ $QV }{ X } || 0),
	 ($calib_stats{ $QV }{ M } || 0) , ($calib_stats{ $QV }{ X } || 0),
      ) if (1);
}

if ( $outfile ) {
  $outfile .= ".bam" if ( $outfile !~ /bam\z/);

  my $reads_recalibrated = 0;

  open (my $bamout, " | $samtools view -o $outfile -Sb - ") || die "Could not open bam stream: $!\n";
#  open (my $bamout, " | echo  ") || die "Could not open bam stream: $!\n";
#  print $bamout `$samtools view -H $bam_file`;
  print $bamout `$samtools view -H $bam_file`;

  open (my $bamin, "$samtools view $bam_file | ") || die "Could not open 'stream': $!\n";
  while(<$bamin>) {

    last if ( $MAX_READS && $MAX_READS <= $reads_recalibrated++);

    my @F = split("\t");
    $F[10] = join("", map{ $calib_matrix{ $_ } || $_} split(//, $F[10]));
    print $bamout join("\t", @F);
  }
  close($bamin);
  close( $bamout );
}

#die Dumper( \%hist );

# 
# 
# 
# Kim Brugger (05 Nov 2010)
sub analyse {
  my ($region) = @_;

  my (@splits, @ref_ids, @reads, @SNPs);

  my $current_pos = undef;

#  print "$samtools view $bam_file $region\n";

  my ($reads_analysed, $bases_analysed) = (0,0);

  open (my $bam, "$samtools view $bam_file $region | ") || die "Could not open 'stream': $!\n";
  while(<$bam>) {

    chomp;
    my @F = split("\t");
    my ($id, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $seq, $qual, @opts) = @F;
    
    # Not a mapped read
    if ($flags & 0x0004 ) {
      next;
    }

    last if ( $MAX_READS && $MAX_READS <= $reads_analysed++);
    last if ( $MAX_BASES && $MAX_BASES < $bases_analysed);

    # pseq: patched sequence, the aligned sequence
    # pqual: the quality values for the aligned bases.
    my ($pseq, $pqual) = patch_alignment($seq, $qual, $cigar);
    my $gseq = substr($fasta, $pos - 1, length( $pseq ));
    $bases_analysed += length($pseq);

    # the cs_splits array needs to be synced with this new pos. Bring forth the 
    # array as many places. If there is a gap, traverse the whole thing and reset
    # the whole thing, so we start from fresh.
    while( @splits && $splits[0]{pos} != $pos ) {
      
#      print "Stepping from $splits[0]{pos} towards $pos\n";
      
      my ($right, $total) = (0, 0);
      map { $total += $splits[0]{ $_ }||0;
	    $right += $splits[0]{ $_ } if ( $splits[0]{ $_ } && $splits[0]{ref} eq $_)}
      ( 'A','C','G','T');
      
      my $ratio = phred_value( $right ,$total - $right);
      my $skip = 0;
      
      
      $skip++	if ( $filter_on_base_mutation_rate && $right != $total && $ratio <= $base_mutation_rate);
      
      if ( ! $skip ) {
	foreach my $base ('A', 'C', 'G', 'T' ) {
	  if ( $splits[0]{ref} eq $base ) {
	    foreach my $QV (keys %{$splits[0]{QV}{$base}} ) {
	    $calib_stats{ $QV }{ $base }{ M } += $splits[0]{QV}{$base}{ $QV };
	    $calib_stats{ $QV }{ M } += $splits[0]{QV}{$base}{ $QV };
	    }
	  }
	  else {
	    foreach my $QV (keys %{$splits[0]{QV}{$base}} ) {
	      $calib_stats{ $QV }{ $base }{ X } += $splits[0]{QV}{$base}{ $QV };
	      $calib_stats{ $QV }{ X } += $splits[0]{QV}{$base}{ $QV };
	    }
	  }
	}
      } 
      
      shift @splits;
      
    }

    my @gseq  = split("", $gseq );
    my @pseq  = split("", $pseq );
    my @pqual = split("", $pqual);
    my $gaps = 0;
    for(my $i = 0; $i<@gseq; $i++ ) {
      
      # there is an insert in the reference, so this number needs to 
      # be subtracted to get the real genome position.
      if ($gseq[ $i ] eq "-") {
	$gaps++;
	next;
      }

      $splits[$i + $gaps]{ ref      } = $gseq[ $i ];
      $splits[$i + $gaps]{ pos      } = $pos  + $i;
      $splits[$i + $gaps]{ $pseq[$i]  }++;
      $splits[$i + $gaps]{ QV }{ $pseq[$i]} {ord($pqual[$i]) - 33  }++;
      next if ($pseq[ $i ] eq "-");
      $splits[$i + $gaps]{ total    }++;
      $hist{ord($pqual[$i]) - 33  }++;
    }
  }

  while( @splits) {
    # the cs_splits array needs to be synced with this new pos. Bring forth the 
    # array as many places. If there is a gap, traverse the whole thing and reset
    # the whole thing, so we start from fresh.
    
    
#    next if (! $splits[0]{total} || $splits[0]{total} == 0);
    
      
    my ($right, $total) = (0, 0);
    map { $total += $splits[0]{ $_ }||0;
	  $right += $splits[0]{ $_ } if ( $splits[0]{ $_ } && $splits[0]{ref} eq $_)}
    ( 'A','C','G','T');
    
    my $ratio = phred_value( $right ,$total - $right);
    my $skip = 0;
    $skip++ if ( $filter_on_base_mutation_rate && $right != $total && $ratio <= $base_mutation_rate);

#    die Dumper( $splits[0]);
      
    if ( ! $skip ) {
      foreach my $base ('A', 'C', 'G', 'T' ) {
	if ( $splits[0]{ref} eq $base ) {
	  foreach my $QV (keys %{$splits[0]{QV}{$base}} ) {
	    $calib_stats{ $QV }{ $base }{ M } += $splits[0]{QV}{$base}{ $QV };
	    $calib_stats{ $QV }{ M } += $splits[0]{QV}{$base}{ $QV };
	  }
	}
	else {
	  foreach my $QV (keys %{$splits[0]{QV}{$base}} ) {
	    $calib_stats{ $QV }{ $base }{ X } += $splits[0]{QV}{$base}{ $QV };
	    $calib_stats{ $QV }{ X } += $splits[0]{QV}{$base}{ $QV };
	  }
	}
      }
    } 
    
    shift @splits;
  }


}


# 
# 
# 
# Kim Brugger (11 Oct 2011)
sub phred_value {
  my ($correct, $wrong) = @_;

#  print "[$wrong] [$correct]\n";

  return -1 if ( !$wrong && ! $correct);
  return 55 if ( ! $wrong && $correct);
  return  0 if ( $wrong && ! $correct);
  
  my $P = $wrong/($wrong+$correct);
  
  return 55 if ( ! $P );
  my $Q = -10 * log10( $P );
  
  return  int ($Q );
}




sub log10 {
  my $n = shift;
  return log($n)/log(10);
}


# 
# 
# 
# Kim Brugger (20 Jul 2009)
sub patch_alignment {
  my ( $seq, $qual, $cigar ) = @_;

  return ($seq, $qual) if ( $cigar !~ /[DIS]/);
  
  my @seq  = split("", $seq );
  my @qual = split("", $qual );


  my (@cigar) = $cigar =~ /(\d*\w)/g;

  my $offset = 0;

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)



  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d+)(\w)/;

    if ( $type eq 'M') {
      $offset += $length;
      next;
    }
    elsif ( $type eq "D") {
      my @dashes = split("", "-"x$length);
      splice(@seq,  $offset, 0, @dashes);
      splice(@qual, $offset, 0, @dashes);
      $offset += $length;
    }
    elsif ( $type eq "I" || $type eq "S" ) {
      splice(@seq,  $offset, $length);
      splice(@qual, $offset, $length);
    }    

  }

  return (join("", @seq), join("",@qual));
}




#
# Read the fasta files and puts entries into a nice array
#
sub readfasta {
  my ($file, $region) = @_;  

  $region =~ s/:\d+-\d+//;

  my $sequence;
  my $header;

  open (my $f, "$samtools faidx $file $region|" ) || die "Could not open $file:$1\n";
  while (<$f>) {
    chomp;
    if (/^\>/) {
      if ($header) { # we have a name and a seq
	return ($header, $sequence);
      }
      $header = $_;
      $header =~ s/^\>//;
    }
    else {$sequence .= $_;}
  }

  return $sequence;
}




# 
# 
# 
# Kim Brugger (13 Jul 2010)
sub find_program {
  my ($program) = @_;
  
  my $username = scalar getpwuid $<;
  
  my @paths = ("/home/$username/bin/",
	       "./",
	       "/usr/local/bin");
  
  foreach my $path ( @paths ) {
    return "$path/$program" if ( -e "$path/$program" );
  }

  my $location = `which $program`;
  chomp( $location);
  
  return $location if ( $location );
  
  die "Could not find '$program'\n";
}


# 
# 
# 
# Kim Brugger (05 Nov 2010)
sub Usage {
  $0 =~ s/.*\///;
  die "USAGE: $0 -b<am file> -R<eference genome (fasta)> -d[ min depth, default=15] -s[ min Split, default=60] -B[uffer, default=100] -M[ set mapq score for offending reads] -U[n set mapped flag for offending reads]\n";

  # tests :::: odd SNP reporting:  10:74879852-74879852
  # large indel: 10:111800742
}


__END__


# 
# 
# 
# Kim Brugger (28 Oct 2010)
sub scrub {
  my ( $read, $SNPs ) = @_;

  return if ( $$read{indel} );

  return if ( ! $$read{singles} && ! $$read{doubles} );
	  
  foreach my $snp ( @$SNPs ) {
    my ($snp_start, $snp_end) = ($$snp[0][0], $$snp[1][0]);
#    next if ( ! $snp_start || ! $snp_end);

    if ( $$read{ pos } - 5 > $snp_end) {
#      print STDERR "Pre-Bailing... $$read{pos} -> $$read{end} vs $snp_start => $snp_end\n";
      $b_post++;
      next;
    }


    if ( $$read{ end } + 5 < $snp_start) {
#      print STDERR "Bailing...\n";
      $b_pre++;
      return;
    }

#    next;
      
    foreach my $single (@{$$read{singles}}) {

      if (($single - 1 + $$read{pos} <= $snp_start &&
	   $single + 1 + $$read{pos} >= $snp_end)
	  || 
	  ($single - 1 + $$read{pos} <= $snp_start &&
	   $single + 1 + $$read{pos} >= $snp_start)
	  ||
	  ($single - 1 + $$read{pos} <= $snp_end &&
	   $single + 1 + $$read{pos} >= $snp_end)) {

	print STDERR "$$read{id} -- $$read{cigar} $$read{flags}\n$$read{a}";
	
	
#	$$read{flags} += 4 if ( $set_unmapped );
#	$$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);
	push @{$$read{sam}}, "ES:Z:1";
	$s++;
	return;
      }
    }  

#    next;

    foreach my $double (@{$$read{doubles}}) {
      
      if (($$double[0] + $$read{pos} <= $snp_start &&
	   $$double[1] + $$read{pos} >= $snp_end)
	  || 
	  ($$double[0] + $$read{pos} <= $snp_start &&
	   $$double[1] + $$read{pos} >= $snp_start)
	  ||
	  ($$double[0] + $$read{pos} <= $snp_end &&
	   $$double[1] + $$read{pos} >= $snp_end)) {
	
#	print STDERR "$$read{id} -- $$read{cigar}\n$$read{a}";
#	print $$read{a};

	$$read{flags} += 4 if ( $set_unmapped );
	$$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);
	push @{$$read{sam}}, "ES:Z:2";
	$d++;
	return;
      }
    }  	  	  
	  
    if ( $$read{pos} == $snp_end || $$read{end}  == $snp_end ) {
	    
#      print STDERR "$$read{id} -- $$read{cigar}\n$$read{a}";
#      print $$read{a};
      $$read{flags} -= 4 if ( $set_unmapped );
      $$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);

      push @{$$read{sam}}, "ES:Z:3";
#      $$read{flags} -= 4;
#      $$read{mapq}   = 2;		
      $e++;
      return;
    }
  }

  return;
}


