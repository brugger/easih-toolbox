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

  my ($res );
  my $errors = "";
  my $warnings = "";
  my $project_name = undef;
  

  open(my $in, $sample_sheet) || die ("Could not open '$sample_sheet': $!\n");
  my @lines;
  while(<$in>) {
    # get rid of all the odd new-lining and quoting of fields.
    $_ =~ s/\r\n/\n/g; 
    $_ =~ s/\n\r/\n/g; 
    $_ =~ s/\r/\n/g; 
    $_ =~ s/\"//g; 
    $_ =~ s/\'//g; 

    push @lines, split("\n",$_);
  }
  close $in;

  if ( $lines[0] =~ /^FCID/) {
    ($res, $errors) = parse_HiSeq_sheet( \@lines );
  }
  elsif ( $lines[0] =~ /^\[Header\]/) {
    ($res, $errors, $project_name) = parse_MiSeq_sheet( \@lines );
  }



  return (undef, $errors, $warnings, $project_name) if ($errors && $errors ne "");

  return ($res, undef, $warnings, $project_name);
}




# 
# 
# 
# Kim Brugger (29 Oct 2012)
sub parse_MiSeq_sheet {
  my ( $lines ) = @_;
  
  my $field_delim = "";
  my ( %res, $project_name );
  my ( $warnings, $errors);

  my $lane = 1;

  while($_ = shift @$lines ) {
    chomp;

    $project_name = $1 if ( $_ =~ /^Project Name\,(.*)/);
    $project_name = $1 if ( ! $project_name && /^Experiment Name,(.*)/);

#    print "$_ || $project_name\n";

    
    # As I dont trust they can export the csv file in the same format each time
    # we will use the first line to identify field and text delimiters.
    if (/^Sample_ID(.)Sample_Name/) {
      $field_delim = $1;

      my ($column_nr, %ids) = (0,());
      map { $ids{ lc($_) } = $column_nr++ } split($field_delim);

      while($_ = shift @$lines ) {
	chomp;
	my @F = split( $field_delim);
	my $barcode = 'default';
	if ( $ids{ index } ) {
	  $barcode = $F[ $ids{ index } ];
	  if ( $ids{ index2 } ) {
	    $barcode .= $F[ $ids{ index2 } ];	    
	  }
	}
	my $sample = $F[ $ids{ lc("Sample_ID") }] || $F[$ids{ lc("Sample_Name") }];
	if ($res{ $lane }{ $barcode } && $res{ $lane }{ $barcode } ne $sample) {
	  if ($barcode eq "default") {
	    $errors .= "Error: Sample without an indexing has already been assigned to '$res{ $lane }{$barcode}' and cannot be assigned to '$sample' as well\n";
	  }
	  else {
	    $errors .=  "Error: Sample with index '$barcode' has already been assigned to '$res{ $lane }{$barcode}' and cannot be assigned to '$sample' as well\n";
	  }
	  next;
	}
	$res{ $lane }{$barcode} = $sample;
      }
    }

  }


  return (\%res, $errors, $project_name);
}

# 
# 
# 
# Kim Brugger (29 Oct 2012)
sub parse_HiSeq_sheet {
  my ( $lines ) = @_;

  my $field_delim = "";
  my %res;
  my ( $warnings, $errors);
  

  while($_ = shift @$lines ) {
    chomp;
    
    # As I dont trust they can export the csv file in the same format each time
    # we will use the first line to identify field and text delimiters.
    if (/^(.{0,1})FCID/) {
      /FCID(.)/;
      $field_delim = $1;
    }      
    else {
      my @F = split($field_delim, $_);
      
      my $fcount = @F; #svvd2 Jan 17 2012
      $warnings .= "Warning: Missing columns - expecting 9 of them...\n", if($fcount < 9); #svvd2 Jan 17 2012
      $warnings .= "Warning: Additional columns found - expecting only 9 of them\n", if($fcount > 9); #svvd2 Jan 17 2012

      
      my (undef, $lane, $sample_id, undef, $index, undef) = @F;

      $index ||= "default";

      $errors .= "Error: Lane number $lane is not a valid lane number - should be between 1 and 8\n", if($lane !~ /^[1-8]$/); #svvd2 18 Jan 2012 


      if ($res{$lane}{$index} && $res{$lane}{$index} ne $sample_id) {
	if ($index eq "default") {
	  $errors .= "Error: Lane $lane without an indexing has already been assigned to '$res{$lane}{$index}' and cannot be assigned to '$sample_id' as well\n";
	}
	else {
	  $errors .=  "Error: Lane $lane with index '$index' has already been assigned to '$res{$lane}{$index}' and cannot be assigned to '$sample_id' as well\n";
	}
	next;
      }

      
      $res{$lane}{$index} = $sample_id;
    }
  }


  return \%res;
  
}




# 
# 
# 
# Kim Brugger (29 Oct 2012)
sub full_MiSeq_sheet {
  my ( $sample_sheet ) = @_;
  
  my ($res );
  my $errors = "";
  my $warnings = "";
  my $project_name = undef;
  
  open(my $in, $sample_sheet) || die ("Could not open '$sample_sheet': $!\n");
  my @lines;
  while(<$in>) {
    # get rid of all the odd new-lining and quoting of fields.
    $_ =~ s/\r\n/\n/g; 
    $_ =~ s/\n\r/\n/g; 
    $_ =~ s/\r/\n/g; 
    $_ =~ s/\"//g; 
    $_ =~ s/\'//g; 

    push @lines, split("\n",$_);
  }
  close $in;


  my $field_delim = "";
  my ( %res  );

  my $lane = 1;

  while($_ = shift @lines ) {
    chomp;

    $project_name = $1 if ( $_ =~ /^Project Name\,(.*)/);
    $project_name = $1 if ( ! $project_name && /^Experiment Name,(.*)/);

#    print "$_ || $project_name\n";

    
    # As I dont trust they can export the csv file in the same format each time
    # we will use the first line to identify field and text delimiters.
    if (/^Sample_ID(.)Sample_Name/) {
      $field_delim = $1;

      my ($column_nr, %ids) = (0,());
      map { $ids{ lc($_) } = $column_nr++ } split($field_delim);

      print Dumper( \%ids );

      while($_ = shift @lines ) {
	chomp;
	my @F = split( $field_delim);
	my ($index1, $index2) = ("","");
	$index1 = $F[ $ids{ index } ]
	    if ( $ids{ index } );
	$index2 .= $F[ $ids{ index2 } ]
	    if ( $ids{ index2 } );

	my $sample_name = $F[ $ids{ lc("Sample_ID") }] || $F[$ids{ lc("Sample_Name") }];

	my $well_name   = $F[ $ids{ lc("sample_well") }];

	$res{ $well_name } = {sample => $sample_name,
			      index1 => $index1,
			      index2 => $index2};

      }
    }

  }


  return (\%res, $errors, $project_name);
}


# 
# 
# 
# Kim Brugger (27 Jul 2011)
sub validate {
  my ($hash, @lanes) = @_;
#  print  Dumper( \$hash );

  my $errors = "";
  
  my %lanes_hash;
  map { $lanes_hash{ $_ } = 2} @lanes;

  foreach my $lane ( @lanes ) {
    
    if (! $$hash{$lane}) {
      $errors .= "Error: no lane information for lane $lane \n";
      next;
    }
   
    if (ref ($$hash{$lane}) eq "HASH") {
      foreach my $bcode (keys %{$$hash{$lane}}) {

	my $sample_name = $$hash{$lane}{$bcode};
#	$sample_name =~ s/-.*//;
	$sample_name =~ s/^(.{7}).*/$1/;
	
	$errors .= "Error: $$hash{$lane}{$bcode} for lane $lane is not an EASIH sample name\n" 
	    if ( !EASIH::Sample::validate_sample( $sample_name ));

	next if ($bcode eq 'default');


	if ( $bcode =~ /^[ACGT]+\z/ ) {
	  die "Error: Both indexed and non-indexed samples in lane $lane\n" if ($$hash{$lane}{'default'});
	  EASIH::Barcodes::barcode_set('illumina');

	  ($bcode, my $valid) = EASIH::Barcodes::validate_barcode( $bcode );
	  $errors .= "Error: '$bcode' in lane $lane is not an valid illumina barcode\n" 
	      if ( ! $valid && length($bcode) != 16 );
	}
	else {
	  $errors .= "Error: '$bcode' in lane $lane is not a valid value, it needs to be a sequence not '$bcode'. \n"; 
	}
      }
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


    $removed{$lane_nr}{ 'default' } = $$sample_sheet{ $lane_nr }{'default'};
    delete( $$sample_sheet{ $lane_nr }{'default'} );
    # this is a mixture of non-index and indexed samples in this lane, remove the non-indexed one as 
    # all the junk will go to this bin.
#    foreach my $barcode ( keys %{$$sample_sheet{ $lane_nr }} ) {
#      next if ( $barcode eq 'default');
#      $removed{$lane_nr}{ $barcode } = $$sample_sheet{ $lane_nr }{ $barcode };
#    }
    
#    $$sample_sheet{$lane_nr} = $$sample_sheet{ $lane_nr }{'default'};
    
  }


  return ($sample_sheet, \%removed);

  exit;
}



1;
