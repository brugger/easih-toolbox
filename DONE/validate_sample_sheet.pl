#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (30 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use Getopt::Std;
use Data::Dumper;

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 1;
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
use EASIH::Mail;
use EASIH::Sample;
use EASIH::Illumina::Sample_sheet;
use EASIH::Barcodes;
EASIH::Barcodes::barcode_set('illumina');

my %opts;
getopts('h', \%opts);

if($opts{h})
{
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

#  $error_message = "";
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


    $error_message .= "No sample sheet present in $indir\n" if ( ! -e $sample_sheet );
    validate_sample_sheet( $sample_sheet );

  }

  SendEmail("sample sheet error in runfolder: $file", $error_message) if ($error_message);
}



# 
# 
# 
# Kim Brugger (30 Jun 2011)
sub validate_sample_sheet {
  my ( $sample_sheet) = @_;

  my ($res, $errors) = EASIH::Illumina::Sample_sheet::readin( $sample_sheet );

  if ( $errors ) {
    $error_message .= $errors;
    return;
  }

  $errors = EASIH::Illumina::Sample_sheet::validate( $res );
  $error_message .= $errors   if ( $errors );
  
}





###########################################################
sub SendEmail {
    my($subject, $message) = @_;

    $subject = "[easih-done] $subject" if ( $subject !~ /easih-dash/);
    EASIH::Mail::send($to, $subject, $message);
}
