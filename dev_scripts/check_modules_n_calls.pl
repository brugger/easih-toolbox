#!/usr/bin/perl -w
# 
#
# Code to test that the modules are sain. When rearranging that part
# of the code, errors often sneak in due to changed module names.
#
# The programs should be able to identify all the programs and modules
# and test that each of the subfunctions are found in the other modules files. 
# All this is order to make sure that renamed modules/functions are kool.
#
# Kim Brugger (Apr 2005), contact: brugger@mermaid.molbio.ku.dk

use strict;
use Data::Dumper;


my $start_dir = "modules/";
my $shared_dir = "shared_modules/";

my %modules = ();
my %functions = ();


&parse_modules($start_dir);
#&parse_modules($shared_dir);
&check_modulecalls("logistics/");


#
# Check all function calls and checks if they are in our list of 
# function calls.
#
# The function also checks that 'require' statements make sense.
#
# Kim Brugger (Apr 2005), contact: brugger@mermaid.molbio.ku.dk

sub check_modulecalls {
  my ($dir) = @_;

  opendir (DIR, $dir) || die "Could not open '$dir': $!\n";
  my @files = readdir DIR;

  foreach my $file (@files) {
    
    next if ($file eq '.' || $file eq '..' || $file =~ /~$/);
    $file = "$dir$file";
    
    # if we are handling a directory, recursion occurs
    if (-d $file) {
      check_modulecalls("$file/");
    }
    else {

      next if ($file =~ /\~/);

      open (FILE, $file) || die "Could not open '$file': $!\n";
      my $sub_function = "";
      while (<FILE>) {
	next if ($file =~ /~$/);
	
	$_ =~ s/#.*//;
	
	if (/^sub (.*?) \{/ || /^sub (.*?)\{/ ) {
	  $sub_function = $1;
	}
	elsif (/(\w+(::\w+)+) \(/ || /(\w+(::\w+)+)\(/) {
	  my $function_call = $1;
    
#          print "$function_call\n";

	  # This looks like a CPAN module, test and see if it is there
#	  next if ($function_call =~ /^[A-Z]/ || $functino);

	  if (! $functions{$function_call}) {
	    print "FUNCTION: $function_call in $file is missing\n";
	  }
	  else {
	    $functions{$function_call}++;
	  }
	}
	elsif (/require (.*)\;/) {
	  my $requirement = $1;

	  next if ($requirement eq "AutoLoader" ||
		   $requirement eq "Exporter");

	  if (! $modules{$requirement}) {
	    if ($sub_function) {
	      print "REQUIRES: $requirement in $file is wrong\n";
	    }
	    else {
	      print "REQUIRES: $requirement in $file, function: $sub_function is wrong\n";
	    }
	  }
	  
	}
      }
    }
  }
}


#
# Traverse the directory tree and check module names/positions
# while we collect a list of all the functions we declare.
#
# Kim Brugger (Apr 2005), contact: brugger@mermaid.molbio.ku.dk
sub parse_modules {
  my ($dir) = @_;

  opendir (DIR, $dir) || die "Could not open '$dir': $!\n";
  my @files = readdir DIR;

  foreach my $file (@files) {
    
    next if ($file eq '.' || $file eq '..' || $file =~ /~$/);
    $file = "$dir$file";
    
    if (-d $file) {
      parse_modules("$file/");
    }
    else {
      open (FILE, $file) || die "Could not open '$file': $!\n";
      my ($package_name) = ("");


      while (<FILE>) {
	my $file2 = $file;
	$file2 =~ s/$shared_dir//;
	$file2 =~ s/$start_dir//;

	next if ($file2 =~ /template.pm/);
	next if ($file2 =~ /\~/);

	if (/^package (.*)\;/) {
	  
	  $package_name = $1;
	  my $hypo_package = "$package_name\.pm";
	  $hypo_package =~ s/::/\//g;
	  
	  if ($hypo_package ne $file2) {
	    print STDERR  "DEFINE: $file is sat to: $hypo_package\n";
	  }
	  else {
	    $modules{$package_name} = 1;
	  }

#	  print "module name $package_name $file\n";
	}
	elsif (/sub (.*?) /) {
	  next if ($1 eq "BEGIN" || $1 eq "END");
	  $functions{"$package_name\::$1"} = 1;
#	  print "$package_name\::$1\n";
	}
      }
    }
  }
}
