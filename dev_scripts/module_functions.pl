#!/usr/bin/perl -w
# 
# Extracts all functions from the various modules.
# 
# 
# Kim Brugger (Apr 2005), contact: brugger@mermaid.molbio.ku.dk

use strict;

my $start_dir = "modules/";

extract_from_dir($start_dir);

sub extract_from_dir {
  my ($dir) = @_;

  opendir (DIR, $dir) || die "Could not open '$dir': $!\n";
  my @files = readdir DIR;

  foreach my $file (@files) {
    
    next if ($file eq '.' || $file eq '..' || $file =~ /~$/);
    $file = "$dir$file";
    
    if (-d $file) {
#      print "Got a dir '$file'\n";
      extract_from_dir("$file/");
    }
    else {
      # check and see if it is a file.
      
      open (FILE, $file) || die "Could not open '$file': $!\n";
      my ($package_name) = ("");

#      print "File: $file\n";

      while (<FILE>) {
	if (/^package (.*)\;/) {
	  $package_name = $1;
	}
	elsif (/sub (.*?) /) {
	  next if ($1 eq "BEGIN" || $1 eq "END");
	  print "$package_name\::$1\n";
#	  print "$package_name\n";
	}
      }
    }
  }
}
