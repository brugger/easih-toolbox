#!/usr/bin/perl 
# 
# splits a fastq file into smaller sub files
# 
# 
# Kim Brugger (22 Apr 2010), contact: kim.brugger@easih.ac.uk


use strict;
use warnings;
use Getopt::Std;


my %opts;
getopts('1:2:ce:n:ho:', \%opts);

my $entries     = $opts{'e'} || 10000000; # 10 mill.
my $first_file  = $opts{'1'} || usage();
my $second_file = $opts{'2'};
my $compress = $opts{c} ? 0 : 1;
my $output_dir  = $opts{o} || "./";
my $outfiles    = $opts{n};
system "mkdir $output_dir" if ( ! -d $output_dir);
usage() if ( $opts{h});

my $file_counter  = 1;
my $entry_counter = 0;
my (@rfh, @wfh);

my $printed = 0;
    
system "mkdir $output_dir" if ( ! -d $output_dir);

if ( $outfiles ) {
  split_to_n_files( );
}
else {
  split_by_count();
}


# 
# 
# 
# Kim Brugger (04 Apr 2011)
sub split_to_n_files {
  
  
  open ($rfh[1], "gzip -dc $first_file |" )  || die "Could not open '$first_file': $!(1)\n"  if ( $first_file  && $first_file =~ /gz/);
  open ($rfh[1], "$first_file" )             || die "Could not open '$first_file': $!(2)\n"  if ( $first_file  && $first_file !~ /gz/);
  open ($rfh[2], "gzip -dc $second_file |" ) || die "Could not open '$second_file': $!(3)\n" if ( $second_file && $second_file =~ /gz/);
  open ($rfh[2], "$second_file" )            || die "Could not open '$second_file': $!(4)\n" if ( $second_file && $second_file !~ /gz/);


  # trim off the absolute path so files does not go somewhere odd.
  $first_file  =~ s/.*\/// if ($first_file);
  $second_file =~ s/.*\/// if ($second_file);
  
  my (@outfiles1, @outfiles2);
  for(my $i=0;$i< $outfiles;$i++) {
    if ( $first_file ) {
      if ( $compress) {
	open ($outfiles1[$i], "| gzip -c  > $output_dir/$first_file.$i.gz")  || die "could not open '$output_dir/$first_file.$file_counter.gz': $!\n";
      }
      else {
	open ($outfiles1[$i], "> $output_dir/$first_file.$i")  || die "could not open '$output_dir/$first_file.$file_counter.gz': $!\n" 
      }
    }

    if ( $second_file ) {
      if ( $compress ) {
	open ($outfiles2[$i], "| gzip -c  > $output_dir/$second_file.$i.gz")  || die "could not open '$output_dir/$second_file.$file_counter.gz': $!\n";
      }
      else {
	open ($outfiles2[$i], "> $output_dir/$second_file.$i")  || die "could not open '$output_dir/$second_file.$file_counter.gz': $!\n";
      }
    }
  }
  
  my ($reads1, $reads2, $reads) = (0,0,0);
  my (@first_counter, @second_counter);
  my @first  = &read1($rfh[1], \$reads1) if ( $first_file ); 
  my @second = &read1($rfh[2], \$reads2) if ( $second_file);
  $reads = $reads1 + $reads2;
  while (@first && @second) {
    
    my $fh = ($reads/2) % $outfiles;
    if ($first[0] eq $second[0]) { # mate pair
      print {$outfiles1[$fh]} $first[1]; 
      print {$outfiles2[$fh]} $second[1];
      
      $first_counter[$fh]++;
      $second_counter[$fh]++;

      @first = &read1($rfh[1], \$reads1); 
      @second = &read1($rfh[2], \$reads2);
      $reads = $reads1 + $reads2;
    } 
    else {
      if ($first[0] le $second[0]) {
	print {$outfiles1[$fh]} $first[1];
	$first_counter[$fh]++;
	@first = &read1($rfh[1], \$reads1);
	$reads = $reads1 + $reads2;
      } 
      else {
	print {$outfiles2[$fh]} $second[1];
	$second_counter[$fh]++;
	@second = &read1($rfh[2]);
	$reads = $reads1 + $reads2;
      }
    }
  }      
  
  if (@first) {
    my $fh = ($reads/2) % $outfiles;
    print {$outfiles1[$fh]} $first[1];
    $first_counter[$fh]++;
    while (@first = &read1($rfh[1], , \$reads1)) {
      $reads = $reads1 + $reads2;
      $fh = ($reads/2) % $outfiles;
      print {$outfiles1[$fh]} $first[1];
      $first_counter[$fh]++;
    }
  }
  
  if (@second) {
    my $fh = $reads % $outfiles;
    print {$outfiles1[$fh]} $second[1];
    $second_counter[$fh]++;
    while (@second = &read1($rfh[1], , \$reads1)) {
      $reads = $reads1 + $reads2;
      $fh = $reads % $outfiles;
      print {$outfiles1[$fh]} $second[1];
      $second_counter[$fh]++;
    }
  }

  use Data::Dumper;
  

  for (my $i=0; $i< @first_counter && $i < @second_counter; $i++) {
    if ( $compress ) {
      print "$output_dir/$first_file.$i.gz\t$output_dir/$second_file.$i.gz\t$first_counter[$i]\t$second_counter[$i]\n";
    }
    else {
      print "$output_dir/$first_file.$i\t$output_dir/$second_file.$i\t$first_counter[$i]\t$second_counter[$i]\n";
    }
  }
  
  print "Total reads:\t$reads\n";
}




# 
# 
# 
# Kim Brugger (04 Apr 2011)
sub split_by_count {
  


  if ( $first_file && $second_file ) {
    
    open ($rfh[1], "gzip -dc $first_file |" )  || die "Could not open '$first_file': $!(1)\n"  if ( $first_file =~ /gz/);
    open ($rfh[1], "$first_file" )             || die "Could not open '$first_file': $!(2)\n"  if ( $first_file !~ /gz/);
    open ($rfh[2], "gzip -dc $second_file |" ) || die "Could not open '$second_file': $!(3)\n" if ( $second_file =~ /gz/);
    open ($rfh[2], "$second_file" )            || die "Could not open '$second_file': $!(4)\n" if ( $second_file !~ /gz/);
    
    # trim off the absolute path so files does not go somewhere odd.
    $first_file  =~ s/.*\///;
    $second_file =~ s/.*\///;
    
    if ( $compress ) {
      open ($wfh[1], "| gzip -c  > $output_dir/$first_file.$file_counter.gz")  || die "could not open '$output_dir/$first_file.$file_counter.gz': $!\n";
      open ($wfh[2], "| gzip -c  > $output_dir/$second_file.$file_counter.gz") || die "could not open '$output_dir/$second_file.$file_counter.gz': $!\n";
      print "$output_dir/$first_file.$file_counter.gz\t$output_dir/$second_file.$file_counter.gz\n";
    }
    else {
      open ($wfh[1], "> $output_dir/$first_file.$file_counter")  || die "could not open '$output_dir/$first_file.$file_counter': $!\n";
      open ($wfh[2], "> $output_dir/$second_file.$file_counter") || die "could not open '$output_dir/$second_file.$file_counter': $!\n";
      print "$output_dir/$first_file.$file_counter\t$output_dir/$second_file.$file_counter\n";
    }  
    $file_counter++;
    
    my (@first, @second);
    @first  = &read1($rfh[1]); 
    @second = &read1($rfh[2]);
    while (@first && @second) {
      
      if ($first[0] eq $second[0]) { # mate pair
	print {$wfh[1]} $first[1]; 
	print {$wfh[2]} $second[1];
	
	@first = &read1($rfh[1]); 
	@second = &read1($rfh[2]);
      } 
      else {
	if ($first[0] le $second[0]) {
	  print {$wfh[1]} $first[1];
	  @first = &read1($rfh[1]);
	} 
	else {
	  print {$wfh[2]} $second[1];
	  @second = &read1($rfh[2]);
	}
      }
      
      $printed++;
      
      if ($printed >= $entries) {
	close($wfh[1]);
	close($wfh[2]);
	
	if ( $compress ) {
	  open ($wfh[1], "| gzip -c > $output_dir/$first_file.$file_counter.gz")  || die "could not open '$output_dir/$first_file.$file_counter.gz': $!\n";
	  open ($wfh[2], "| gzip -c > $output_dir/$second_file.$file_counter.gz") || die "could not open '$output_dir/$second_file.$file_counter.gz': $!\n";
	  print "$output_dir/$first_file.$file_counter.gz\t$output_dir/$second_file.$file_counter.gz\n";
	}
	else {
	  open ($wfh[1], "> $output_dir/$first_file.$file_counter")  || die "could not open '$output_dir/$first_file.$file_counter': $!\n";
	  open ($wfh[2], "> $output_dir/$second_file.$file_counter") || die "could not open '$output_dir/$second_file.$file_counter': $!\n";
	  print "$output_dir/$first_file.$file_counter\t$output_dir/$second_file.$file_counter\n";
	}  
	
	$file_counter++;
	$printed = 0;
      }
    }
    
    
    if (@first) {
      print {$wfh[1]} $first[1];
      while (@first = &read1($rfh[1])) {
	print {$wfh[1]} $first[1];
      }
    }
    
    if (@second) {
      print {$wfh[2]} $second[1];
      while (@second = &read1($rfh[2])) {
	print {$wfh[2]} $second[1];
      }
    }
  }
  else {
    
    open ($rfh[1], "gzip -dc $first_file |" )  || die "Could not open '$first_file': $!\n"  if ( $first_file =~ /gz/);
    open ($rfh[1], "$first_file" )             || die "Could not open '$first_file': $!\n"  if ( $first_file !~ /gz/);
    
    $first_file  =~ s/.*\///;
    
    if ( $compress ) {
      open ($wfh[1], "| gzip -c > $output_dir/$first_file.$file_counter.gz")  || die "could not open '$output_dir/$first_file.$file_counter.gz': $!\n";
      print "$output_dir/$first_file.$file_counter.gz\n";
    }
    else {
      open ($wfh[1], "> $output_dir/$first_file.$file_counter")  || die "could not open '$output_dir/$first_file.$file_counter': $!\n";
      print "$output_dir/$first_file.$file_counter\n";
    }  
    $file_counter++;
    
    while ( my @first  = &read1($rfh[1])) {
      
      print {$wfh[1]} $first[1];
      $printed++;
      
      
      if ($printed >= $entries) {
	close($wfh[1]);
	
	if ( $compress ) {
	  open ($wfh[1], "| gzip -c > $output_dir/$first_file.$file_counter.gz")  || die "could not open '$output_dir/$first_file.$file_counter.gz': $!\n";
	  print "$output_dir/$first_file.$file_counter.gz\n";
	}
	else {
	  open ($wfh[1], "> $output_dir/$first_file.$file_counter")  || die "could not open '$output_dir/$first_file.$file_counter': $!\n";
	  print "$output_dir/$first_file.$file_counter\n";
	}  
	
	$file_counter++;
	$printed = 0;
      }
    }
  }
}


# 
# 
# 
# Kim Brugger (06 Aug 2010)
sub read1 {
  my ( $read_handle, $counter) = @_;

  my $key = <$read_handle>;
  return () if ( ! $key );
  my $seq = <$read_handle>;
  my $strand = <$read_handle>;
  my $qual   = <$read_handle>;
  
  die "fastq file is either not a fq-file or contains extra lines.\n"  
      if ( $key !~ /^\@/ || $strand !~ /[+-]/ );

  my $entry = $key . $seq . $strand. $qual;
  chomp($key);
  #chomp off /1|/2 id from the solid names
  $key =~ s/\/[1|2|3]//;
  $$counter++;
  return $key ? ($key, $entry) : ();
}


sub usage {

  $0 =~ s/.*\///;
  print "USAGE: $0 --entries <number, def=10mill>  -o[ut dir, ./ is default] -1 fq-file -2 fq-file\n";
  exit;
}
