package EASIH::Illumina::Sample_sheet;
# 
# Illumina sample sheet handling module
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use EASIH::Barcodes;

BEGIN { 
  EASIH::Barcodes::barcode_set('illumina');
}

# 
# Examines a sample_sheet_hash to determine if it is a barcoded run.
# 
# Kim Brugger (13 Sep 2011)
sub indexed_run {
  my ( $sample_hash ) = @_;
  
  print Dumper( $sample_hash );
  
  foreach my $lane ( keys %{$sample_hash})  {
    if ( ref($$sample_hash{ $lane }) eq "HASH" &&
	 !$$sample_hash{ $lane }{'default'} ) {
      return 1;
    }
  }

  return 0;
}



# 
# 
# 
# Kim Brugger (27 Jul 2011)
sub readin {
  my ( $sample_sheet ) = @_;

  my (%res );
  my $errors = "";
  
  my $text_delim = "";
  my $field_delim = "";

  open(my $in, $sample_sheet) || die ("Could not open '$sample_sheet': $!\n");
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

      $index ||= "default";

      $lane      =~ s/^$text_delim(.*)$text_delim\z/$1/;
      $sample_id =~ s/^$text_delim(.*)$text_delim\z/$1/;
      $index     =~ s/^$text_delim(.*)$text_delim\z/$1/;

      if ($res{$lane}{$index} && $res{$lane}{$index} ne $sample_id) {
	if ($index eq "") {
	  $errors .= "Lane $lane without an indexing has already been assigned to '$res{$lane}{$index}' and cannot be assigned to '$sample_id' as well\n";
	}
	else {
	  $errors .=  "Lane $lane with index '$index' has already been assigned to '$res{$lane}{$index}' and cannot be assigned to '$sample_id' as well\n";
	}
	next;
      }


      $res{$lane}{$index} = $sample_id;
    }
  }

  foreach my $lane ( keys %res ) {
    if (keys %{$res{$lane}} == 1 && $res{$lane}{'default'}) {
      $res{$lane} = $res{$lane}{'default'};
    }
  }

  return (undef, $errors) if ($errors ne "");

  return (\%res, undef);
}



# 
# 
# 
# Kim Brugger (27 Jul 2011)
sub validate {
  my ($hash, $limited_lanes) = @_;
#  print  Dumper( \$hash );

  my $errors = "";
  

  for ( my $lane =1; $lane <=8;$lane++) {
    
    if (! $$hash{$lane} && !$limited_lanes) {
      $errors .= "no lane information for lane $lane \n";
      next;
    }

    
   
    if (ref ($$hash{$lane}) eq "HASH") {
      foreach my $bcode (keys %{$$hash{$lane}}) {
	
	$errors .= "$$hash{$lane}{$bcode}} for lane $lane is not an EASIH sample name\n" 
	    if ( !EASIH::Sample::validate_name($$hash{$lane}{$bcode}));

	next if ($bcode eq 'default');


	if ( $bcode =~ /^[ACGT]+\z/ ) {
	  die "Both indexed and non-indexed samples in lane $lane\n" if ($$hash{$lane}{'default'});
	  EASIH::Barcodes::barcode_set('illumina');

	  ($bcode, my $valid) = EASIH::Barcodes::validate_barcode( $bcode );
	  $errors .= "'$bcode' in lane $lane is not an valid illumina barcode\n" if ( ! $valid );
	}
	elsif($bcode =~ /^(\w+?)_(\w+)([F|R]{1,2})/) {
	  my $barcode_set = $1;
	  my $tag = $2;
	  my $directions = $3;

	  die "Only supports ill9 EASIH barcodes, not $barcode_set\n" if ( $barcode_set ne "ill9");
	  die "Only supports m13 cloning site, not '$tag'\n" if ( $tag ne "m13");


	}
	else {
	  $errors .= "'$bcode' in lane $lane is not a valid value. It should be either be an illumina barcode or an EASIH barcode group\n"; 
	}
      }
    }
    else {
	$errors .= "$$hash{$lane} for lane $lane is not an  EASIH sample name\n" if ( !EASIH::Sample::validate_name($$hash{$lane}));
    }
  }


  return undef if ($errors eq "");
  return $errors;
}



# 
# 
# 
# Kim Brugger (20 Sep 2011)
sub remove_easih_barcodes {
  my ( $sample_sheet ) =  @_;


  my %removed;


  foreach my $lane_nr ( keys %$sample_sheet ) {

    next if ( ref ($$sample_sheet{ $lane_nr }) ne "HASH");
    next if ( ref ($$sample_sheet{ $lane_nr }) eq "HASH" && ! $$sample_sheet{ $lane_nr }{'default'});

    foreach my $barcode ( keys %{$$sample_sheet{ $lane_nr }} ) {
      next if ( $barcode eq 'default');
      $removed{$lane_nr}{ $barcode } = $$sample_sheet{ $lane_nr }{ $barcode };
    }
    
    $$sample_sheet{$lane_nr} = $$sample_sheet{ $lane_nr }{'default'};
    
  }


  return ($sample_sheet, \%removed);

  exit;
}



1;
