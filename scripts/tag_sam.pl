#!/usr/bin/perl 
# 
# Add or replace a readgroup/sample/library/platform information to a samfile.
# 
# 
# Kim Brugger (22 Jul 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use File::Temp qw/ tempfile /;


my %opts;
getopts('i:o:Or:s:l:p:c:a:A:R:', \%opts);
usage() if ( $opts{h});

my $infile    = $opts{i};
my $outfile   = $opts{o};
my $readgroup = $opts{r} || usage();
$readgroup    =~ s/.*\///;
my $reference = $opts{R} || "";
$reference =~ s/.*\///;
$reference =~ s/^(.*?)\..*/$1/;

my $sample    = $opts{s} || $readgroup;
$sample =~ s/.*\///;
my $library   = $opts{l} || $readgroup;
$library =~ s/.*\///;
my $platform  = $opts{p} || usage();
my $center    = $opts{c} || "EASIH";

my $aligner   = $opts{a};
my $a_line    = $opts{A};

my $replace   = $opts{O};
$infile = $replace if ( ! $infile && $replace );
open (*STDIN, $infile) || die "Could not open '$infile': $!\n" if ( $infile );

my ($tmp_fh, $tmp_file);
if ( $replace) {
  ($tmp_fh, $tmp_file) =  tempfile( DIR => './');
  *STDOUT = *$tmp_fh;
}
elsif ( $outfile ) {
  open (*STDOUT, " > $outfile " ) || die "Could not open '$outfile': $!\n" ;
}


my ($written_readgroup, $written_cmd_flags) = (0, 0);

while(<STDIN>) {
  
  if ( /^\@RG/) {
    my ( $same_readgroup, $same_sample, $same_library, $same_platform, $seq_center) = (0, 0, 0, 0);
    foreach my $field ( split("\t", $_)) {
      next if ($field =~ /^\@/);
      my ($key, $value) = split(":", $field);
      $same_readgroup++ if ( $field eq 'ID' && $value eq $readgroup);
      $same_sample++    if ( $field eq 'SM' && $value eq $sample);
      $same_library++   if ( $field eq 'LB' && $value eq $library);
      $same_platform++  if ( $field eq 'PL' && $value eq $platform);
      $seq_center++     if ( $field eq 'CN' && $value eq $center);
    }

#    $written_readgroup++ if ( $same_readgroup && $same_sample && $same_library && $same_platform );
    if ( $same_readgroup && $same_sample && $same_library && $same_platform && $seq_center ) {
      $written_readgroup++;
    }
    else {
      next;
    }
  }

  if ( /^\@PG/) {
    my ( $same_aligner, $same_param) = (0, 0);
    foreach my $field ( split("\t", $_)) {
      next if ($field =~ /^\@/);
      my ($key, $value) = split(":", $field);
      $same_aligner++ if ( $field eq 'ID' && $value eq $aligner);
      $same_param++    if ( $field eq 'CL' && $value eq $a_line);
    }

    if ( $same_aligner && $same_param ) {
      $written_cmd_flags++;
    }
    else {
      next;
    }
  }

  if ($reference && /^\@SQ/ && !/\tAN:/) {
    chomp;
    $_ .= "\tAN:$reference\n";
  }
    
    

  if ( ! /^\@/ ) {

    if ( $aligner && ! $written_cmd_flags ) {
      my @fields = ("\@PG", "ID:$aligner");
      push @fields, "CL:$a_line" if ( $a_line );
      print join("\t", @fields) . "\n";
      $written_cmd_flags++;
    }

    if (! $written_readgroup ) {
      print join("\t", "\@RG", "ID:$readgroup","SM:$sample","LB:$library","PL:$platform", "CN:$center") . "\n";
      $written_readgroup++;
    }

  
    if ( (/\tRG:Z:(\w+)\t/ || /\tRG:Z:(\w+)\Z/) && 
	 (/\tSM:Z:(\w+)\t/ || /\tSM:Z:(\w+)\Z/)) {
      
      s/(.*\tRG:Z:)(.*?)(\t.*)/$1$readgroup$3/;
      s/(.*\tSM:Z:)(.*?)(\t.*)/$1$sample$3/;
      
      s/(.*\tRG:Z:)(.*?)\Z/$1$readgroup/;
      s/(.*\tSM:Z:)(.*?)\Z/$1$sample/;
    }
    elsif ( /\tRG:Z:(\w+)\t/ || /\tRG:Z:(\w+)\Z/ ) {
      chomp($_);
      s/(.*\tRG:Z:)(.*?)(\t.*)/$1$readgroup$3/;
      s/(.*\tRG:Z:)(.*?)\Z/$1$readgroup/;
      $_ .= "\tSM:$sample\n";
    }
    elsif ( /\tSM:Z:(\w+)\t/ || /\tSM:Z:(\w+)\Z/ ) {
      chomp($_);
      s/(.*\tSM:Z:)(.*?)(\t.*)/$1$sample$3/;
      s/(.*\tSM:Z:)(.*?)\Z/$1$sample/;
      $_ .= "\tRG:$readgroup\n";
    }
    else {
      chomp($_);
      $_ .= "\tRG:Z:$readgroup\tSM:Z:$sample\n";
    }
  }

  print;
  
}


if ( $replace) {
  close($tmp_fh);
  system "mv $tmp_file $infile";
}




# 
# 
# 
# Kim Brugger (22 Jul 2010)
sub usage {

  $0 =~ s/.*\///;
  print "$0 adds/replaces the readgroup/library/sample/platform/center tags in a sam file/stream\n";
  print "$0 -i[nfile (or stdin] -o[utfile (or stdout)] -r[eadgroup] -s[ample] -l[ibrary] -p[latform] -c[enter] -a[ligner] -A[ligner param] -O[verwrite the original file] -R[eference]\n";

  exit 1;
}

