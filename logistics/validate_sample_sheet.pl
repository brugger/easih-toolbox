#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Jun 2011), contact: kim.brugger@easih.ac.uk

use lib "/home/kb468/easih-toolbox/modules/";
use strict;
use Getopt::Std;
use Data::Dumper;

use EASIH::Mail;
use EASIH::Logistics;

my %opts;
getopts('h', \%opts);

if($opts{h})
{
    $0 =~ s/.*\///;
    print STDERR "\n\nDescription: Checks for sample sheet in the most recent illumina run folder, and validates that it is correct\n";
    print STDERR "\n\nDescription: if not correct, will send an notification email\n";
    exit -1;
}

my $error_message = "";


#my $to = 'sri.deevi@easih.ac.uk,kim.brugger@easih.ac.uk'; #global
my $to = 'kim.brugger@easih.ac.uk'; #global
#my $to = 'bics@easih.ac.uk,lab@easih.ac.uk'; #global
my @input_dirs = ("/seqs/illumina2/");

foreach my $dir ( @input_dirs ) {

  opendir(DIR, "$dir");
  my @files = grep(!/^\.|\.log$/, sort readdir(DIR));
  closedir(DIR);

  my $file = pop @files;

  $error_message = "";

  if ( -e "$dir$file/Data/Intensities/BaseCalls") {

    chomp(my $checkstring = `egrep -c RunStart_Netcopy $dir$file/Events.log`);
    
    next if ( ! $checkstring);
    my $indir = "$dir$file/Data/Intensities/BaseCalls";
    my $sample_sheet = "";
    $sample_sheet = "$indir/sample_sheet.csv" if (!$sample_sheet && -e "$indir/sample_sheet.csv");
    $sample_sheet = "$indir/Sample_sheet.csv" if (!$sample_sheet && -e "$indir/Sample_sheet.csv");
    $sample_sheet = "$indir/sample_Sheet.csv" if (!$sample_sheet && -e "$indir/sample_Sheet.csv");
    $sample_sheet = "$indir/Sample_Sheet.csv" if (!$sample_sheet && -e "$indir/Sample_Sheet.csv");


#    print "Examining $dir$file $sample_sheet\n";


    fail("No sample sheet present in $indir\n") if ( ! -e $sample_sheet );
    
    validate_sample_sheet( $sample_sheet );
  }


  SendEmail("sample sheet error in runfolder: $file", $error_message);
  
}



# 
# 
# 
# Kim Brugger (30 Jun 2011)
sub validate_sample_sheet {
  my ( $sample_sheet) = @_;

  my (%res );
  
  my $text_delim = "";
  my $field_delim = "";

  open(my $in, $sample_sheet) || fail("Could not open '$sample_sheet': $!\n", "BASECALL2FQ_PATH_ERROR");
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

      $index ||= "";

      $lane      =~ s/^$text_delim(.*)$text_delim\z/$1/;
      $sample_id =~ s/^$text_delim(.*)$text_delim\z/$1/;
      $index     =~ s/^$text_delim(.*)$text_delim\z/$1/;

      fail( "Index should be a base sequence, not '$index' for lane $lane\n")  if ( $index && $index !~ /^[ACGT]\z/i);

      if ( $index ) {
	fail( "Lane $lane with index '$index' has already been assigned to '$res{$lane}{$index}' and cannot be assigned to '$sample_id' as well\n") 
	    if ($res{$lane}{$index} && !$opts{$lane} && !$opts{'a'});

	$res{$lane}{$index} = $sample_id;
      }
      else {
	fail( "Lane $lane has already been assigned to '$res{$lane}' and cannot be assigned to '$sample_id' as well\n") 
	    if ($res{$lane}  && !$opts{$lane} && !$opts{'a'});

	$res{$lane} = $sample_id;
      }
    }
  }

  my %basenames;
  for ( my $lane =1; $lane <=8;$lane++) {
    
    fail( "no lane information for lane $lane \n")
	if (! $res{$lane});
    
  
    if (ref ($res{$lane}) eq "HASH") {
      foreach my $bcode (keys %{$res{$lane}}) {
	fail("$res{$lane}{$bcode}} for lane $lane is not an EASIH sample name\n") if ( !EASIH::Logistics::validate_sample_name($res{$lane}{$bcode}));
      }
    }
    else {
	fail("$res{$lane} for lane $lane is not an  EASIH sample name\n") if ( !EASIH::Logistics::validate_sample_name($res{$lane}));
    }
  }
  
}




# 
# 
# 
# Kim Brugger (30 Jun 2011)
sub fail {
  my ($message) = @_;

  $error_message .= "$message";
  
}


###########################################################
sub SendEmail {
    my($subject, $message) = @_;
 
    $subject = "[easih-dash] $subject" if ( $subject !~ /easih-dash/);
    EASIH::Mail::send($to, $subject, $message);
}
