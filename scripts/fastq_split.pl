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
getopts('1:2:e:nho:', \%opts);

my $entries     = $opts{'e'} || 10000000; # 10 mill.
my $first_file  = $opts{'1'} || usage();
my $second_file = $opts{'2'};
my $compress = $opts{n} ? 0 : 1;
my $output_dir  = $opts{o} || "./";
system "mkdir $output_dir" if ( ! -d $output_dir);
usage() if ( $opts{h});

my $file_counter  = 1;
my $entry_counter = 0;
my (@rfh, @wfh);

my $printed = 0;

if ( $first_file && $second_file ) {

  open ($rfh[1], "gzip -dc $first_file |" )  || die "Could not open '$first_file': $!\n"  if ( $first_file =~ /gz/);
  open ($rfh[1], "$first_file" )             || die "Could not open '$first_file': $!\n"  if ( $first_file !~ /gz/);
  open ($rfh[2], "gzip -dc $second_file |" ) || die "Could not open '$second_file': $!\n" if ( $second_file =~ /gz/);
  open ($rfh[2], "$second_file" )            || die "Could not open '$second_file': $!\n" if ( $second_file !~ /gz/);

  # trim off the absolute path so files does not go somewhere odd.
  $first_file  =~ s/.*\///;
  print "--> $first_file\n";
  $second_file =~ s/.*\///;

  system "mkdir $output_dir" if ( ! -d $output_dir);
  
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


# 
# 
# 
# Kim Brugger (06 Aug 2010)
sub read1 {
  my ( $read_handle, $match) = @_;

  my $key = <$read_handle>;
  return () if ( ! $key );
  my $entry = $key . <$read_handle> . <$read_handle>. <$read_handle>;
  chomp($key);
  #chomp off /1|/2 id from the solid names
  $key =~ s/\/[1|2]//;
  return $key ? ($key, $entry) : ();
}


sub usage {

  $0 =~ s/.*\///;
  print "USAGE: $0 --entries <number, def=10mill>  -o[ut dir, ./ is default] -1 fq-file -2 fq-file\n";
  exit;
}
