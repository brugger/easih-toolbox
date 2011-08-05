package EASIH::Fastq;
# 
# Generic fastq IO
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


# 
# 
# Kim Brugger (03 Aug 2011)
sub open {
  my ($filename) = @_;

  my $fh;

  if ( $filename =~ /gz/) {
    CORE::open ( $fh, "gunzip -c $filename | ") || die "Could not open '$filename': $!\n";
  }
  else {
    CORE::open ( $fh, "$filename") || die "Could not open '$filename': $!\n";
  }

  return $fh;
  
}



# 
# 
# 
# Kim Brugger (03 Aug 2011)
sub next {
  my ($fh) = @_;

  my $name = <$fh>;
  my $seq = <$fh>;  
  my $strand  =  <$fh>;
  my $qual = <$fh>;

  return 0  if (! $seq);

  return ($name, $seq, $strand, $qual);
  
}




1;



