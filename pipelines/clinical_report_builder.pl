#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 May 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Spreadsheet::WriteExcel;

use lib '/software/packages/easih-toolbox/modules';
use EASIH::Git;

my $PIPELINE_VERSION = '1.4';


my $variations;
my %columns;
my $variation_file = shift;
my $easih_id = $variation_file;
#$easih_id =~ s/^\.\.\///;
$easih_id =~ s/\..*//;
my $excel_file .= "$easih_id.xls";
my $leeway = 30;

# Create a new Excel workbook
my $workbook = Spreadsheet::WriteExcel->new( $excel_file );
my $bold     = $workbook->add_format(bold => 1);

my %added_worksheets;
my %worksheet_offset;

my %common_polys;
$common_polys{BRCA1}{'c.1-134'} = "T>C";
$common_polys{BRCA1}{'c.301+26'} = "C>A";
$common_polys{BRCA1}{'c.301+27'} = "A>C";
$common_polys{BRCA1}{'c.442-34'} = "T>C";
$common_polys{BRCA1}{'c.548-58'} = "delT";
$common_polys{BRCA1}{'c.1067'} = "A>G";
$common_polys{BRCA1}{'c.2077'} ="G>A";
$common_polys{BRCA1}{'c.2082'} ="C>T";
$common_polys{BRCA1}{'c.2311'} ="T>C";
$common_polys{BRCA1}{'c.2612'} ="C>T";
$common_polys{BRCA1}{'c.3113'} ="A>G";
$common_polys{BRCA1}{'c.3119'} ="G>A";
$common_polys{BRCA1}{'c.3548'} ="A>G";
$common_polys{BRCA1}{'c.3708'} ="T>G";
$common_polys{BRCA1}{'c.4308'} ="T>C";
$common_polys{BRCA1}{'c.4837'} ="A>G";
$common_polys{BRCA2}{'c.1-26'} = "G>A";
$common_polys{BRCA2}{'c.865'} = "A>C";
$common_polys{BRCA2}{'c.1114'} = "C>A";
$common_polys{BRCA2}{'c.1365'} = "A>G";
$common_polys{BRCA2}{'c.2229'} = "T>C";
$common_polys{BRCA2}{'c.2971'} = "A>G";
$common_polys{BRCA2}{'c.3396'} = "A>G";
$common_polys{BRCA2}{'c.3807'} = "T>C";
$common_polys{BRCA2}{'c.7242'} = "A>G";
$common_polys{BRCA2}{'c.7806-14'} = "C>T";
$common_polys{BRCA2}{'c.8567'} = "A>C";
$common_polys{BRCA2}{'c.8851'} = "G>A";
$common_polys{BRCA2}{'c.9649-19'} = "G>A";
$common_polys{BRCA2}{'c.9257-16'} = "T>C";
$common_polys{BRCA2}{'c.9746'} = "C>T";

my @col_order = ('exon', 'Position', 'Nucleotide pos', 'Change', 
		 'AA change', 'Score', 'Depth', 'AAF', #'FP',
		 'Common poly', 'dbsnp', 'HGMD', 'Comment', 'Repeat','Check1', 'Check2');

$added_worksheets{ 'QC' } = $workbook->add_worksheet('QC');

$added_worksheets{ 'QC' }->write(0, 0, 'EASIH ID');
$added_worksheets{ 'QC' }->write(0, 1, $easih_id, $bold);
$added_worksheets{ 'QC' }->write(1, 0, 'GM No');
$added_worksheets{ 'QC' }->write(0, 3, 'Pipeline version:');
$added_worksheets{ 'QC' }->write(0, 4, "$PIPELINE_VERSION-".EASIH::Git::version(), $bold);

$added_worksheets{ 'QC' }->write(4, 0, 'Name', $bold); 
$added_worksheets{ 'QC' }->write(4, 1, 'Min depth', $bold);     
$added_worksheets{ 'QC' }->write(4, 2, 'Max depth', $bold);     
$added_worksheets{ 'QC' }->write(4, 3, 'Mean depth', $bold);    
$added_worksheets{ 'QC' }->write(4, 4, 'Lows', $bold);    
$added_worksheets{ 'QC' }->write(4, 5, 'Missing', $bold);
$worksheet_offset{ 'QC' } = 5;

#exit;
my @header;
my %vars;

open(my $var_fh, $variation_file) || die "Could not open '$variation_file': $!\n";
while(<$var_fh>) {
  # Ignore all lines with comments.
  if (/#\s+capture: (.*)/) {
#    print "$1\n";
    chomp($1);
    my @F = split("\t", $1);
    for(my $i=0;$i<@F;$i++) {
      $F[$i] =~ s/'//g;
      $added_worksheets{ 'QC' }->write($worksheet_offset{ 'QC' }, $i, $F[$i]);
    }
      $worksheet_offset{ 'QC' }++;
    next;
  }
  elsif( /\#/ ) {
    next;
  }

  chomp;

  # This is the line that describes what is in what column, pull this
  # information out and index it and we will use it later for
  # accessing the correct data from each line.  this makes the code
  # maore flexible if the columns change...
  if ( /^Position/) {
    my @column_names = split("\t");
    @header = split("\t");
    my $column_index  = 0;
    foreach my $column_name ( @column_names) {
      $columns{ $column_name } = $column_index++;
    }
    next;
  }
  
  my @values = split("\t");

  my $gene    = $values[$columns{ 'gene' }];
  my $exon    = $values[$columns{ 'exon' }];

  next if ($exon ne "INTRON" && $exon !~ /$gene/);

  my $effect  = $values[$columns{ 'Effect' }];;

#  print "$effect\n";

  next if ( $effect =~ /UTR/ || $effect =~ /STREAM/);

  if ( ! $added_worksheets{ $gene }) {

    $added_worksheets{ $gene } = $workbook->add_worksheet($gene);
    # add some std cells to all work sheets.
    $added_worksheets{ $gene }->write(4, 0, 'EASIH ID');
    $added_worksheets{ $gene }->write(4, 1, $easih_id, $bold);
    $added_worksheets{ $gene }->write(5, 0, 'GM No');
    $added_worksheets{ $gene }->write(4, 6, 'Gene');
    $added_worksheets{ $gene }->write(4, 7, $gene );
    $added_worksheets{ $gene }->write(5, 6, 'Transcript');
    $added_worksheets{ $gene }->write(5, 7, $values[$columns{ 'transcript' }] );
    
    $worksheet_offset{ $gene } = 10;
    for(my $i=0;$i<@col_order;$i++) {
      $added_worksheets{ $gene }->write($worksheet_offset{ $gene }, $i, 
					ucfirst($col_order[$i]), $bold);
    }
    $worksheet_offset{ $gene }++;

    $added_worksheets{ "$gene full" } = $workbook->add_worksheet("$gene full");
    # add some std cells to all work sheets.
    $added_worksheets{ "$gene full" }->write(4, 0, 'EASIH ID');
    $added_worksheets{ "$gene full" }->write(4, 1, $easih_id, $bold);
    $added_worksheets{ "$gene full" }->write(5, 0, 'GM No');
    $added_worksheets{ "$gene full" }->write(4, 6, 'Gene');
    $added_worksheets{ "$gene full" }->write(4, 7, $gene );
    $added_worksheets{ "$gene full" }->write(5, 6, 'Transcript');
    $added_worksheets{ "$gene full" }->write(5, 7, $values[$columns{ 'transcript' }] );
    
    $worksheet_offset{ "$gene full" } = 10;
    for(my $i=0;$i<@col_order;$i++) {

      $added_worksheets{ "$gene full" }->write($worksheet_offset{ "$gene full" }, $i, 
					ucfirst($col_order[$i]), $bold);
    }
    $worksheet_offset{ "$gene full" }++;

  }
  
  my $cpos = $values[$columns{ 'Nucleotide pos' }];
  
  $vars{ $gene }{ $cpos } = $_;
}


#print Dumper (\%vars );

foreach my $gene (keys %vars ) {
  foreach my $cpos ( sort { pure_cpos($a) <=> pure_cpos($b)  || pure_cpos_offset($a) <=> pure_cpos_offset($b)} keys %{$vars{ $gene }} ) {
    
    next if ($cpos =~ /\*/);

    my @values = split("\t", $vars{$gene}{$cpos});


    my $pure_cpos   = pure_cpos($cpos);
    my $cpos_offset = pure_cpos_offset($cpos);

#    print "$cpos $pure_cpos $cpos_offset <= $leeway $cpos\n";

    for(my $i=0;$i<@col_order;$i++) {
      my $col_name = $col_order[$i];
      my $col_nr = $columns{ $col_name };
      my $col_value = "";
      $col_value = $values[$col_nr] if ( defined $col_nr && $values[$col_nr] );
      $col_value = "INTRON" if ( $col_value eq "" && $col_name eq "exon");
      $col_value = "Yes" if ( $col_name eq "Common poly" && $common_polys{$gene}{$cpos});


      if ($col_name eq "exon" && $col_value =~ /exon/i) {
	$col_value =~ s/ INTRON//;
      }


      if ($col_name eq "dbsnp" && $col_value ) {
	my $rs_number = $col_value;
	$rs_number =~ s/rs//;
	$rs_number = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=$rs_number";

	$added_worksheets{ "$gene full" }->write($worksheet_offset{ "$gene full" }, $i,  $rs_number, $col_value);
	$added_worksheets{ $gene }->write($worksheet_offset{ $gene }, $i,  $rs_number, $col_value)
	    if ($cpos_offset <= $leeway);
      }
      else {
	$added_worksheets{ "$gene full" }->write($worksheet_offset{ "$gene full" }, $i,  $col_value);
	$added_worksheets{ $gene }->write($worksheet_offset{ $gene }, $i,  $col_value)
	    if ($cpos_offset <= $leeway);
      }
    }
    
    
    $worksheet_offset{ "$gene full" }++;
    $worksheet_offset{ $gene }++ if ($cpos_offset <= $leeway);
    
  }
}


for my $sheet (keys %added_worksheets) {

  next if $sheet eq 'QC';

  $worksheet_offset{$sheet} += 5;

  my $i = 0;
  my @headings = ("Exon", "PCR W/S", "Seq W/S", "ABI Run", "Result", "PCR W/S", "Seq W/S", "ABI Run", "Result", "Check1", "Check2");
  for my $head (@headings) {
    $added_worksheets{$sheet}->write($worksheet_offset{$sheet}, $i, $head, $bold);
    $i++;
  }

}


# 
# 
# 
# Kim Brugger (29 Jun 2012)
sub pure_cpos {
  my $pos = shift;

  
  $pos =~ s/del.*\z//;
  $pos =~ s/ins.*\z//;

  $pos =~ s/\s//g;

  $pos =~ s/c\.//;
  $pos =~ s/^\*//;
  $pos =~ s/\+\d+\z//;
  $pos =~ s/\-\d+\z//;
  $pos =~ s/\s+//g;
  return 0 if ( $pos eq "");
  return $pos;
		
  
}


# 
# 
# 
# Kim Brugger (29 Jun 2012)
sub pure_cpos_offset {
  my $pos = shift;

  my $track = 0;
  $track = 1 if ( $pos =~ /68-638/);
  


  $pos =~ s/del.*\z//;
  $pos =~ s/ins.*\z//;
  $pos =~ s/\s//g;

  return 0 if ( $pos !~ /\+\d+\z/ && $pos !~ /\-\d+\z/);



  $pos =~ s/c\.//;
  $pos =~ s/^\*//;
  $pos =~ s/^.*\+//;
  $pos =~ s/^.*\-//;
  $pos =~ s/\s+//g;


  return 0 if ( $pos eq "");
  return $pos;
		
  
}


