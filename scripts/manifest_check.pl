#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (16 Mar 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::MD5;


if ( ! -e 'MANIFEST') {
  print STDERR "No MANIFEST file present in current directory cannot preform check\n";
  exit -1;
}

open (my $mani, "MANIFEST") || die "Could not open 'MANIFEST': $!\n";

my %manifest_files;
map {chomp;    
     next if ($_ =~ /^\z/);
     $manifest_files{$_}++ } <$mani>;
close( $mani );

my %local_files;
open(my $f, "find ./ |") || die "Could not open find pipeline: $!\n";
while (<$f>) {
  chomp;
  s/\.\///;
  next if ($_ =~ /^\z/);
  if ( ! $manifest_files{ $_ }  && -f $_ ) {
    print STDERR "The file '$_' is present in the directory but don't appear in the MANIFEST\n";
    exit -1;
  }

  if ( ! $manifest_files{ $_ }  && -d $_ ) {
    print STDERR "The directory '$_' is present in the directory but don't appear in the MANIFEST\n";
    exit -1;
  }
  
  $local_files{$_}++;
}
close($f);


foreach my $file ( keys %manifest_files ) {
  if (! $local_files{ $file } ) {
    print "'$file' present in the MANIFEST is missing\n";
    exit -1;
  }
  if ( $file =~ /md5/) {
    if (! EASIH::MD5::validate_file( $file ) ) {
      print "md5 for '$file' is wrong\n";
      exit -1;
    }
    else {
#      print "md5 for '$file' is correct\n";
    }
  }
  
}


print "All files in the MANIFEST are present and all the md5 sums are correct\n";



# 
# 
# 
# Kim Brugger (16 Mar 2012)
sub local_files {
  
}
