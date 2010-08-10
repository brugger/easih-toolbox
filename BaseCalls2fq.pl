#!/usr/bin/perl 
# 
# Transforms a illumina basecalls folder into 8 or 16 fq files depending on 
# wether is was a single or a paired run. 
# 
# Kim Brugger (10 Aug 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

my %opts;
getopts("a:1:2:3:4:5:6:7:8:hs:d:o:", \%opts);

my $lane1 = $opts{'1'};
my $lane2 = $opts{'2'};
my $lane3 = $opts{'3'};
my $lane4 = $opts{'4'};
my $lane5 = $opts{'5'};
my $lane6 = $opts{'6'};
my $lane7 = $opts{'7'};
my $lane8 = $opts{'8'};

$lane1 = $opts{a} if ($opts{a});
$lane2 = $opts{a} if ($opts{a});
$lane3 = $opts{a} if ($opts{a});
$lane4 = $opts{a} if ($opts{a});
$lane5 = $opts{a} if ($opts{a});
$lane6 = $opts{a} if ($opts{a});
$lane7 = $opts{a} if ($opts{a});
$lane8 = $opts{a} if ($opts{a});

my @lanes = (undef, $lane1, $lane2, $lane3, $lane4, $lane5, $lane6, $lane7, $lane8);
my $indir   = $opts{'d'} || "./";
my $outdir  = $opts{'o'} || "./";
system "mkdir $outdir" if ( ! -d $outdir);

#parse_SampleSheet ( $opts{s} );

validate_lane_names();

for(my $i = 1; $i<=8; $i++) {

  my $lane_name = $lanes[ $i ];

  my @files = glob("$indir/s_$i\_1_*_qseq.txt");
  if ( @files ) {
    open (my $out, "| gzip -c > $out_dir/$lane_name.1.fq.gz") || die "Could not open '$out_dir/$lane_name.1.fq.gz': $!\n";
    my ($count_in, $count_out) = analyse_files($out, @files);
    printf ("lane $i.1\t$lane_name\t$count_in\t$count_out (%.2f %%)\t%.2f avg clusters per tile\n", $count_out*100/$count_in, $count_out/120) ;
  }

  @files = glob("$indir/s_$i\_2_*_qseq.txt");
  if ( @files ) {
    open (my $out, "| gzip -c > $out_dir/$lane_name.2.fq.gz") || die "Could not open '$out_dir/$lane_name.2.fq.gz': $!\n";
    my ($count_in, $count_out) = analyse_files($out, @files);
    printf ("lane $i.2\t$lane_name\t$count_in\t$count_out (%.2f %%)\t%.2f avg clusters per tile\n", $count_out*100/$count_in, $count_out/120) ;
  }
  last;
}



# 
# 
# 
# Kim Brugger (10 Aug 2010)
sub analyse_files {
  my ($fhw, @files) = @_;

  my ($count_in, $count_out) = (0,0);

  foreach my $file ( @files) {
    open (my $in, "$file") || die "Could not open '$file': $!\n";
    
    while (my $line = <$in>) {
      
      chomp $line;
      $count_in++;

      my ($instr, $run_id, $lane, $tile, $x, $y, $index, $read, $bases, $q_line, $filter)
	  = split /\t/, $line;

      #Did not pass the chastity filter
      next if ( $filter == 0);

      $bases =~ tr/./N/;           # turn dots into Ns
      $q_line =~ tr/!-\175/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-\136/;
      
      if ($index ne '0') {
        print $fhw "\@${instr}_$run_id:$lane:$tile:$x:$y\#$index/$read\n";
      } else {
        print $fhw "\@${instr}_$run_id:$lane:$tile:$x:$y/$read\n";
      }
      print $fhw $bases, "\n";
      print $fhw "+\n";
      print $fhw $q_line, "\n";
      $count_out++;
      
    }
  }

  return ($count_in, $count_out);
}



# 
# 
# 
# Kim Brugger (10 Aug 2010)
sub validate_lane_names {

  my (%seen, $error);
  my @missing_names;
  for( my $i=1; $i< @lanes; $i++) {
    my $name = $lanes[ $i ];
    if ( ! defined $name || $name eq "") {
      push @missing_names, $i;
      next;
    }
    push @{$seen{$name}}, $i;
  }

  die "Missing names for the following lanes: @missing_names\n" if ( @missing_names);

  
  
  foreach my $name ( keys %seen) {
    
    if ( @{$seen{$name}} > 1)  {
      my @new_names;
      
      map {$lanes[ $_ ] .= "_$_"; push @new_names,$lanes[ $_ ]} @{$seen{$name}};
      print "$name is being used multiple times. Renaming $name to @new_names \n";
    }
  }
}

