package EASIH::Sample_sheet;
# 
# EASIH sample sheet handling module
# format of sheet: Lane\tSample (w/ path & version)\tMID/barcode
# 
# Kim Brugger (25 May 2012), contact: kim.brugger@easih.ac.uk


use strict;
use warnings;
use Data::Dumper;

use EASIH::Barcodes;

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
  my %sample; #svvd2 Jan 17 2012
  
  # convert lines as they could come from either Linux, mac or Windows.
  open(my $in, $sample_sheet) || die ("Could not open '$sample_sheet': $!\n");
  my @lines;
  while(<$in>) {
    $_ =~ s/\r\n/\n/g; 
    $_ =~ s/\n\r/\n/g; 
    $_ =~ s/\r/\n/g; 
    $_ =~ s/\"//g; 
    push @lines, split("\n",$_);
  }
  close $in;

  while($_ = shift @lines ) {
    chomp;
    $_ =~ s/\"//g; 
    
    # As I dont trust they can export the csv file in the same format each time
    # we will use the first line to identify field and text delimiters.
    if (/^Lane/) {
      next;
    }      
    elsif( /^(\d+)\t(.+?)\t(\w*)\z/) {
      my ($lane, $sample, $mid) = ($1, $2, $3);

      $mid = 'default' if (  !$mid || $mid eq "");

      if ($res{$lane}{$mid} && $res{$lane}{$mid} ne $sample) {
	if ($mid eq "") {
	  $errors .= "Error: Lane $lane without an indexing has already been assigned to '$res{$lane}{$mid}' and cannot be assigned to '$sample' as well\n";
	}
	else {
	  $errors .=  "Error: Lane $lane with index '$mid' has already been assigned to '$res{$lane}{$mid}' and cannot be assigned to '$sample' as well\n";
	}
	next;
      }

      $res{$lane}{$mid} = $sample;
    }
    else {
      $errors .= "Not a header or illegally formatted line: '$_'\n";
      next;
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





1;
