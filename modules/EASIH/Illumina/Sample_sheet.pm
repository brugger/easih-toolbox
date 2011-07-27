package EASIH::Illumina::Sample_sheet;
# 
# Illumina sample sheet handling module
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;




# 
# 
# 
# Kim Brugger (27 Jul 2011)
sub readin {
  my ( $sample_sheet) = @_;

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

  return (undef, $errors) if ($errors ne "");

  return (\%res, undef);
}



# 
# 
# 
# Kim Brugger (27 Jul 2011)
sub validate {
  my ($hash) = @_;
#  print  Dumper( \$hash );

  my $errors = "";
  
  foreach my $lane ( keys %$hash ) {
    if (keys %{$$hash{$lane}} == 1 && $$hash{$lane}{'default'}) {
      $$hash{$lane} = $$hash{$lane}{'default'};
    }
  }

  for ( my $lane =1; $lane <=8;$lane++) {
    
    if (! $$hash{$lane}) {
      $errors .= "no lane information for lane $lane \n";
      next;
    }
   
    if (ref ($$hash{$lane}) eq "HASH") {
      foreach my $bcode (keys %{$$hash{$lane}}) {
	
	$errors .= "$$hash{$lane}{$bcode}} for lane $lane is not an EASIH sample name\n" 
	    if ( !EASIH::Sample::validate_name($$hash{$lane}{$bcode}));

	next if ($bcode eq 'default');


	if ( $bcode =~ /^[ACGT]+\z/ ) {
	  ($bcode, my $valid) = EASIH::Barcodes::validate_barcode( $bcode );
	  $errors .= "$bcode  in lane $lane is not an valid illumina barcode\n" if ( ! $valid );
	}
	elsif($bcode =~ /^(\w+?)_(\w+)([F|R]{1,2})/) {
	  my $barcode_set = $1;
	  my $tag = $2;
	  my $directions = $3;
	  
	  print "$barcode_set barcodes using the $tag tag for direction(s): $directions\n";
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



1;



