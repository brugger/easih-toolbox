package EASIH::Misc;
#
# Misc functions that does not fit anywhere else, but multiple scripts rely on. 
# 
# 
# 
# Kim Brugger (13 Jul 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use File::Temp;



# 
# 
# 
# Kim Brugger (13 Jul 2010)
sub find_program {
  my ($program, @paths) = @_;

  my $username = scalar getpwuid $<;

  my $location = `which $program`;
  chomp( $location);
  
  return $location if ( $location );

  push @paths, "/home/$username/bin/";
  push @paths, "/home/$username/easih-toolbox/scripts/";
  
  foreach my $path ( @paths ) {
    return "$path/$program" if ( -e "$path/$program" );
  }

  return undef;
}





# 
# 
# 
# Kim Brugger (03 Feb 2011)
sub tmp_dir_file {
  my ( $dir ) = @_;

  my $tmp_dir = tmp_dir( $dir );
  
  my $tmp_file = tmp_file( "." );
  
  return( $tmp_dir, $tmp_file);
  
}


sub tmp_dir {
  my ($dir) = @_;
  
  my ($tmp_fh, $tmp_dir) = File::Temp::tempfile(DIR => "$dir" ) if ( $dir );
  ($tmp_fh, $tmp_dir) = File::Temp::tempfile( ) if ( ! $dir );
  close ($tmp_fh);

  system "rm -f $tmp_dir";
  mkdir "$tmp_dir";
  return "$tmp_dir";
}


sub tmp_file {
  my ($tmp_dir) = @_;

  $tmp_dir ||= '/tmp/';
  
  my ($tmp_fh, $tmp_file) = File::Temp::tempfile(DIR => "$tmp_dir" ) if ( $tmp_dir );
  close ($tmp_fh);
  system "rm -f $tmp_file";

  $tmp_file =~ s/^\.\///;

  return "$tmp_file";
}



# 
# 
# 
# Kim Brugger (16 Feb 2012)
sub random_string {
  my ($length) = @_;

  $length ||= 10;

  my @charset = (('A'..'Z'), ('a'..'z'),(0..9));
  my $string = "";
  while( $length--) {
    $string .= $charset[int(rand(@charset))];
  }
  
  return $string;
}




1;
