package EASIH::Barcodes;
# 
# Generic demultiplexing module. How it works? Well lets find out
# 
# 
# Kim Brugger (25 Jul 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use POSIX qw(ceil floor);

use EASIH::Barcodes::ill9;


my $bc_module   = undef;
my $bc_length   = 0;
my %barcode_hash;
my $bcodes;
my $strict_tags = 1;
my $trim_bc     = 1;
my $correct_bc  = 0;

my $m13fwd         = "TGTAAAACGACGGCCAGT";
my $m13fwd_loose   = "TGTAA.*?CAGT";
my $m13rev         = "CAGGAAACAGCTATGACC";
my $m13rev_loose   = "CAGGA.*?GACC";


# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub error_correct_barcodes {
  my ($correct) = @_;

  $correct_bc = $correct if (defined $correct);

  return $correct_bc;
}


# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub barcode_set {
  my ($bcs) = @_;

  $bcs =~ s/EASIH::Barcodes:{0,2}//;
  $bc_module  = "EASIH::Barcodes::$bcs" if ($bcs);


  $bc_length    = $bc_module->bc_length();
  %barcode_hash = $bc_module->barcodes();
  $bcodes     = [sort keys %barcode_hash];

  return $bc_module;
}


# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub strict_tags {

  $strict_tags = 1;
}


# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub loose_tags {

  $strict_tags = 0;
}



# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub trim_off_bc {
  my ($trim) = @_;
  $trim_bc = $trim;
}



# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub decode_m13f {
  my ( $seq, $qual ) = @_;
  
  if ( ! $strict_tags ) {
    return generic_decode($m13fwd_loose, $seq, $qual);
  }
  else {
    return generic_decode($m13fwd, $seq, $qual);
  }

}


# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub decode_m13r {
  my ( $seq, $qual ) = @_;
  
  if ( ! $strict_tags ) {
    return generic_decode($m13rev_loose, $seq, $qual);
  }
  else {
    return generic_decode($m13rev, $seq, $qual);
  }

}




# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub generic_decode {
  my ($tag_seq,  $seq, $qual ) = @_;

  my ($bc, $res_seq, $res_qual) = (undef, $seq, $qual);

  if ( $seq =~ m/^(.{$bc_length})($tag_seq)(.*)/) {
    
    my ($bcode, $tag, $seq, $found) = ( $1, $2, $3, 0);

    ($bcode, $found) = validate_barcode($bcode);
    
    return (undef, $found, $seq, $qual) if (!$found);

    if ( $trim_bc ) {
      $res_seq = $seq;
      $res_qual = substr($qual, length($bcode)+length($tag));
    }
    
    return ($barcode_hash{$bcode}, $found, $seq, $qual);
  }

  return (undef, -1, $seq, $qual);
  
}



# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub validate_barcode {
  my ($barcode) = @_;

  return ($barcode, 1) if ( $barcode_hash{ $barcode });
  return ($barcode, 0) if ( ! $correct_bc );


  my $START = 0;
  my $END   = 1;
  
#  return 0 if ( $pos != 57421888);


  # set the start and end of the array and find the 
  # the middle of the array
  my ( $left, $right ) = (0, int(@$bcodes)-1);
  my $middle = floor(($right - $left)/2);

  while (1) {

    verbose("MIDDLE $middle ( $left, $right)\n", 1);
    verbose(" $barcode <=> $$bcodes[ $left ] $$bcodes[ $middle ] $$bcodes[ $right ] \n", 1);
    
    if ( one_seq_error($barcode, $$bcodes[ $left ])) {
#      print "Corrected: $barcode == $$bcodes[ $left ]\n";
      return ($$bcodes[ $left ], 2);
    }
    elsif ( one_seq_error($barcode, $$bcodes[ $right ])) {
#      print "Corrected: $barcode == $$bcodes[ $right ]\n";
      return ($$bcodes[ $right ], 2);
    }

    # The new block is to the left of the middle.
    elsif ( $barcode lt $$bcodes[ $middle ] ) {
      $right = $middle;
      $middle = $left + floor(($right - $left)/2);
      verbose("L");
      last if ( $right <= $left || $middle == $left || $middle == $right);
    }
    # The new block is to the right of the middle.
    elsif ($barcode gt $$bcodes[ $middle ] ) {
      $left = $middle;
      $middle = $left + floor(($right - $left)/2);
      verbose("R");
      last if ( $right <= $left || $middle == $left || $middle == $right);
    }
    elsif ( $barcode eq $$bcodes[ $middle ]) {
#      print "Validated: $barcode == $$bcodes[ $middle ]\n";
      return ($barcode, 1);
    }
    elsif ( one_seq_error($barcode, $$bcodes[ $middle ])) {
#      print "Corrected: $barcode == $$bcodes[ $middle ]\n";
      return ($$bcodes[ $middle ], 2);
    }
  }

#  print "$barcode\n";

  return ($barcode, 0);
}


# 
# 
# 
# Kim Brugger (17 May 2011)
sub one_seq_error {
  my ($seq1, $seq2) = @_;

  my @seq1 = split('', $seq1);
  my @seq2 = split('', $seq2);

  my $diffs = 0;
  for(my $i = 0; $i < @seq1; $i++) {
    $diffs++ if ( $seq1[$i] ne $seq2[$i]);
  }

  return 1 if ( $diffs < 2);
  return 0;
}


sub verbose {
  return;
  my ($message, $level) = @_;
  $level  ||= 1;
  $message =~ s/\n+\Z//g;
  print STDERR "MESS " . ":"x$level . " $message\n";
#  sleep 1;
}





# 
# 
# 
# Kim Brugger (25 Jul 2011)
sub BEGIN {
  barcode_set('ill9');
  
}

1;



