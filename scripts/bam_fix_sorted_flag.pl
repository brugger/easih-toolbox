#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (14 Jul 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use File::Temp qw/ tempfile /;

my $bamfile = shift;

if (! $bamfile ) {
  $0 =~ s/.*\///;
  print "USAGE: $0 bam-file. This will set the SO flag to coordinate if the SO flag is set in the original header.\n";
  exit 1;
}

my ($header_tmp_fh, $header_tmp_file) =  tempfile( DIR => './');
my ($bam_tmp_fh, $bam_tmp_file) =  tempfile( DIR => './');
close($bam_tmp_fh);

my $samtools = find_program('samtools');

open (my $bh, "$samtools view -H $bamfile | " ) || die "Could not open samtools stream ($bamfile): $!\n";
while (<$bh>) {
  s/SO:queryname/SO:coordinate/;
  print $header_tmp_fh $_;
}

close($header_tmp_fh);


system "samtools reheader $header_tmp_file $bamfile > $bam_tmp_file.bam";
system "mv -f $bam_tmp_file.bam $bamfile";

system "rm $header_tmp_file";




# 
# 
# 
# Kim Brugger (13 Jul 2010)
sub find_program {
  my ($program) = @_;


  my @paths = ("/home/easih/bin/",
	       "/home/kb468/bin/",
	       "/usr/local/bin");
  
  foreach my $path ( @paths ) {
    
    return "$path/$program" if ( -e "$path/$program" );
  }

  my $location = `which $program`;
  chomp( $location);
  
  return $location if ( $location );

  return undef;
}
