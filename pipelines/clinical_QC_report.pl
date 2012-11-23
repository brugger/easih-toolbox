#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (08 Nov 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use Spreadsheet::WriteExcel;


# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 0;
  if ( $DYNAMIC_LIB_PATHS ) {
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules" && ! grep /^$path\/modules/, @INC);
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    push @INC, '/home/kb468/easih-toolbox/modules/';
  }
}

use EASIH;
use EASIH::Illumina::Sample_sheet;

my %opts;
#getopts("a:A:1:2:3:4:5:6:7:8:hs:Si:o:lhmb:L:", \%opts);
getopts("hB:L:bfSds:", \%opts);

my $sample_sheet = $opts{'s'};

print "$sample_sheet\n";

my ($sample_names, $errors, $warnings, $project_name ) = EASIH::Illumina::Sample_sheet::full_MiSeq_sheet( $sample_sheet );

print Dumper($sample_names);

my $cwd = `pwd`;
chomp( $cwd );

$cwd =~ /(CP\d+)/;

die "Could not resolve the pool id: $cwd\n" if ( ! $1 );

my $excel_file = "$1.xls";

my $workbook = Spreadsheet::WriteExcel->new( $excel_file );
my $bold     = $workbook->add_format(bold => 1);
my $bold_center = $workbook->add_format(bold => 1, valign  => 'vcenter', 
					align   => 'center', border =>1, 
					bg_color     => 'grey',);

my $red_cell  = $workbook->add_format(bg_color     => 'red',);
my $green_cell  = $workbook->add_format(bg_color     => 'green',);

my $QC_tab   = $workbook->add_worksheet('QC');

$QC_tab->merge_range('A2:A7', 'A', $bold_center);
$QC_tab->merge_range('A8:A13', 'B', $bold_center);
$QC_tab->merge_range('A14:A19', 'C', $bold_center);
$QC_tab->merge_range('A20:A25', 'D', $bold_center);
$QC_tab->merge_range('A26:A31', 'E', $bold_center);
$QC_tab->merge_range('A32:A37', 'F', $bold_center);
$QC_tab->merge_range('A38:A43', 'G', $bold_center);
$QC_tab->merge_range('A44:A49', 'H', $bold_center);

$QC_tab->merge_range('B1:C1', '1', $bold_center);
$QC_tab->merge_range('D1:E1', '2', $bold_center);
$QC_tab->merge_range('F1:G1', '3', $bold_center);
$QC_tab->merge_range('H1:I1', '4', $bold_center);
$QC_tab->merge_range('J1:K1', '5', $bold_center);
$QC_tab->merge_range('L1:M1', '6', $bold_center);
$QC_tab->merge_range('N1:O1', '7', $bold_center);
$QC_tab->merge_range('P1:Q1', '8', $bold_center);
$QC_tab->merge_range('R1:S1', '9', $bold_center);
$QC_tab->merge_range('T1:U1', '10', $bold_center);
$QC_tab->merge_range('V1:W1', '11', $bold_center);
$QC_tab->merge_range('X1:Y1', '12', $bold_center);

my $format_top_left    = $workbook->add_format(top => 1, left=>1);
my $format_left        = $workbook->add_format(left=>1);
my $format_left_bottom = $workbook->add_format(left=>1, bottom=>1);
my $format_bottom = $workbook->add_format(bottom=>1);


my ( $ROW, $COLUMN) = (0,1);

my %well_starts = ( "A01" => [1,1],
		    "A02" => [1,3],
		    "A03" => [1,5],
		    "A04" => [1,7],
		    "A05" => [1,9],
		    "A06" => [1,11],
		    "A07" => [1,13],
		    "A08" => [1,15],
		    "A09" => [1,17],
		    "A10" => [1,19],
		    "A11" => [1,21],
		    "A12" => [1,23],

		    "B01" => [7,1],
		    "B02" => [7,3],
		    "B03" => [7,5],
		    "B04" => [7,7],
		    "B05" => [7,9],
		    "B06" => [7,11],
		    "B07" => [7,13],
		    "B08" => [7,15],
		    "B09" => [7,17],
		    "B10" => [7,19],
		    "B11" => [7,21],
		    "B12" => [7,23],


		    "C01" => [13,1],
		    "C02" => [13,3],
		    "C03" => [13,5],
		    "C04" => [13,7],
		    "C05" => [13,9],
		    "C06" => [13,11],
		    "C07" => [13,13],
		    "C08" => [13,15],
		    "C09" => [13,17],
		    "C10" => [13,19],
		    "C11" => [13,21],
		    "C12" => [13,23],

		    "D01" => [19,1],
		    "D02" => [19,3],
		    "D03" => [19,5],
		    "D04" => [19,7],
		    "D05" => [19,9],
		    "D06" => [19,11],
		    "D07" => [19,13],
		    "D08" => [19,15],
		    "D09" => [19,17],
		    "D10" => [19,19],
		    "D11" => [19,21],
		    "D12" => [19,23],

		    "E01" => [25,1],
		    "E02" => [25,3],
		    "E03" => [25,5],
		    "E04" => [25,7],
		    "E05" => [25,9],
		    "E06" => [25,11],
		    "E07" => [25,13],
		    "E08" => [25,15],
		    "E09" => [25,17],
		    "E10" => [25,19],
		    "E11" => [25,21],
		    "E12" => [25,23],

		    "F01" => [31,1],
		    "F02" => [31,3],
		    "F03" => [31,5],
		    "F04" => [31,7],
		    "F05" => [31,9],
		    "F06" => [31,11],
		    "F07" => [31,13],
		    "F08" => [31,15],
		    "F09" => [31,17],
		    "F10" => [31,19],
		    "F11" => [31,21],
		    "F12" => [31,23],

		    "G01" => [37,1],
		    "G02" => [37,3],
		    "G03" => [37,5],
		    "G04" => [37,7],
		    "G05" => [37,9],
		    "G06" => [37,11],
		    "G07" => [37,13],
		    "G08" => [37,15],
		    "G09" => [37,17],
		    "G10" => [37,19],
		    "G11" => [37,21],
		    "G12" => [37,23],

		    "H01" => [43,1],
		    "H02" => [43,3],
		    "H03" => [43,5],
		    "H04" => [43,7],
		    "H05" => [43,9],
		    "H06" => [43,11],
		    "H07" => [43,13],
		    "H08" => [43,15],
		    "H09" => [43,17],
		    "H10" => [43,19],
		    "H11" => [43,21],
		    "H12" => [43,23],);


my %nextera_xt = ('TAAGGCGA' => 'N701',
		  'CGTACTAG' => 'N702',
		  'AGGCAGAA' => 'N703', 
		  'TCCTGAGC' => 'N704', 
		  'GGACTCCT' => 'N705', 
		  'TAGGCATG' => 'N706', 
		  'CTCTCTAC' => 'N707', 
		  'TAGATCGC' => 'N501',
		  'CTCTCTAT' => 'N502',
		  'TATCCTCT' => 'N503',
		  'AGAGTAGA' => 'N504',
		  'GTAAGGAG' => 'N505',
		  'ACTGCATA' => 'N506',
		  'AAGGAGTA' => 'N507',
		  'CTAAGCCT' => 'N508',
		  'CAGAGAGG' => 'N708',
		  'GCTACGCT' => 'N709',
		  'CGAGGCTG' => 'N710',
		  'AAGAGGCA' => 'N711',
		  'GTAGAGGA' => 'N712',);


foreach my $well ( keys %well_starts ) {
  my ($start_row, $start_col) = @{$well_starts{ $well }};
  
  my $sample_name = $$sample_names{ $well }{ 'sample' } || "";
  my $index1      = $$sample_names{ $well }{ 'index1' } || "";
  my $index2      = $$sample_names{ $well }{ 'index2' } || "";

  
  print "$well -- $start_row - $start_col\n";
  
  $QC_tab->write($start_row, $start_col, 'Sample name:', $format_top_left);
  $QC_tab->write($start_row, $start_col+1, $sample_name);


  $QC_tab->write($start_row+1, $start_col, 'index1:', $format_left);
  $QC_tab->write($start_row+1, $start_col+1, $nextera_xt{$index1} || $index1);

  $QC_tab->write($start_row+2, $start_col, 'index2:', $format_left);
  $QC_tab->write($start_row+2, $start_col+1, $nextera_xt{$index2} || $index2);

  my $flagstats = flagstats("$sample_name.bam.flagstat");


  $QC_tab->write($start_row+3, $start_col, 'nr reads:', $format_left);
  $QC_tab->write($start_row+3, $start_col+1, $$flagstats{total} || "");

  $QC_tab->write($start_row+4, $start_col, 'reads mapped:', $format_left);
  $QC_tab->write($start_row+4, $start_col+1, $$flagstats{mapped}-$$flagstats{dups} || "");

  $QC_tab->write($start_row+5, $start_col, 'report done:', $format_left_bottom);

  if ( $sample_name  eq "") {
     $QC_tab->write($start_row+5, $start_col+1, "", $format_bottom);
  }
  elsif ( -e "$sample_name.xls" ) {
     $QC_tab->write($start_row+5, $start_col+1, "Yes", $green_cell, $format_bottom);
  }
  else {
    print "Failed finding $sample_name.xls\n";
     $QC_tab->write($start_row+5, $start_col+1, "No", $red_cell, $format_bottom);
  }
    
  


}  
  
  
  
  



# 
# 
# 
# Kim Brugger (24 Sep 2012)
sub flagstats {
  my ($file ) = @_;

  return undef if ( !-e $file);

  my %res;

  open (my $f, $file) || die "Could not open '$file': $!\n";
  while(<$f>) {
    chomp;

    if ( /(\d+) \+ 0 in total/ ) {
      $res{ total } = $1;
    }
    elsif ( /(\d+) \+ 0 duplicates/ ) {
      $res{ dups } = $1;
    }
    elsif ( /(\d+) \+ 0 mapped/ ) {
      $res{ mapped } = $1;
    }

  }
  close( $f);

  return \%res;
}
