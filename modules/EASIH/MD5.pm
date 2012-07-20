package EASIH::MD5;
# 
# simple md5 module, needs to be extended.
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use EASIH;


# 
# 
# 
# Kim Brugger (13 Jan 2012)
sub create_file {
  my ($infile, $outfile) = @_;
  
  $outfile = "$infile.md5" if ( !$outfile );


  print "md5sum $infile\n";
  
  open( my $p, "md5sum $infile | ") || die "Could not create md5sum pipe: $!\n";
  my $md5 = <$p>;
  close( $p );
  ($md5, $infile) = split(/\s+/, $md5);
  print "$md5 $infile\n";
  $infile =~ s/.*\///;
  open( my $o, "> $outfile ") || die "Could not write to '$outfile': $!\n";
  print $o "$md5  $infile\n";
  close( $o );

  return $md5;
}



# 
# 
# 
# Kim Brugger (13 Jan 2012)
sub validate_file {
  my ($md5file, $infile) = @_;

  open (my $i, $md5file) || die "Could not open '$md5file': $!\n";
  my ($md5sum, $file) = split(/\s+/, <$i>);
  close( $i );

  $infile = $file if (! $infile );

  # set the path to the position of the md5 file.
  $infile = "$1$infile" if ($md5file =~ /(.*\/)/);

  if (validate_sum($md5sum, $infile)) {
    return 1;
    print "Good md5\n";
  }
  else {
    return 0;
    print "Bad md5\n";
  }
}


# 
# 
# 
# Kim Brugger (13 Jan 2012)
sub validate_sum {
  my ($md5sum, $infile) = @_;
  
  open( my $p, "md5sum $infile | ") || die "Could not create md5sum pipe: $!\n";
  my $md5 = <$p>;
  close( $p );
  ($md5, undef) = split(/\s+/, $md5);
  
  return 1 if ( $md5sum eq $md5);

  return 0;
}



1;



