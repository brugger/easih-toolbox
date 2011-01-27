package EASIH::QC;
#
# QC functions that I can hopefully share between projects...
# 
# 
# 
# Kim Brugger (26 Jan 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $sample_size = 0;


# 
# If doing random sampling the number of MB to look at. If -1 or 0 the whole dataset is used 
# 
# Kim Brugger (26 Jan 2011)
sub sample_size {
  my ($size) = @_;
  $sample_size = $size*1048576 if (defined  $size );

  return $sample_size;
}


# 
# Does simple QC on a fastq file without mapping it... If doing a PE run, call the function twice...
# 
# Kim Brugger (03 Aug 2010)
sub csfastQC  {
  my ( $infile1, $res ) = @_;
  
  my %res;
  my $read = 0;
  my $file1;
  if ( $infile1 =~ /gz/) {
    open ( $file1, "gunzip -c $infile1 | ") || die "Could not open '$infile1': $!\n";
  }
  else {
    open ( $file1, "$infile1") || die "Could not open '$infile1': $!\n";
  }
  while (<$file1>) { 
    next if (/#/ || /\>/);
    chomp; 
    my $sequence = substr($_, 1);
    $read += length( $sequence );
    last if ( $sample_size > 0 && $read > $sample_size );
#    print "$read > $sample_size\n";
#    print "$sequence\n";
    $res = analyse( $sequence, undef, $res);
  }

  return $res;
}

# Kim Brugger (03 Aug 2010)
sub qualQC  {
  my ( $infile1, $res ) = @_;
  
  my %res;
  my $read = 0;
  my $file1;
  if ( $infile1 =~ /gz/) {
    open ( $file1, "gunzip -c $infile1 | ") || die "Could not open '$infile1': $!\n";
  }
  else {
    open ( $file1, "$infile1") || die "Could not open '$infile1': $!\n";
  }
  while (<$file1>) { 
    next if (/#/ || /\>/);
    chomp; 

    s/^(\d+)\s*//;
    s/-1/0/g;
    my $qual;
    map { $qual .= chr($_+33)} split(/\s+/);
    
    $read += length( $qual );


    last if ( $sample_size > 0 && $read > $sample_size );
#    print "$read > $sample_size\n";
#    print "$sequence, $strand, $quality\n";
    $res = analyse( undef, $qual, $res);
  }

  return $res;
}
  

# 
# Does simple QC on a fastq file without mapping it... If doing a PE run, call the function twice...
# 
# Kim Brugger (03 Aug 2010)
sub fastQC {
  my ( $infile1) = @_;

  my %res;
  my $read = 0;
  my $file1;
  if ( $infile1 =~ /gz/) {
    open ( $file1, "gunzip -c $infile1 | ") || die "Could not open '$infile1': $!\n";
  }
  else {
    open ( $file1, "$infile1") || die "Could not open '$infile1': $!\n";
  }
  while (<$file1>) { 
#    print "-- $_\n";
    my $sequence = <$file1>;
    my $strand  = <$file1>;
    my $quality  = <$file1>;
    chomp( $sequence); 
    chomp( $strand );
    chomp( $quality);
    $read += length( $sequence );
    last if ( $sample_size > 0 && $read > $sample_size );
#    print "$read > $sample_size\n";
#    print "$sequence, $strand, $quality\n";
    analyse( $sequence, $quality, \%res);
  }

  return \%res;
}


# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub bamQC {
  my ( $infile) = @_;

  my (%res1, %res2);
  my ($read1, $read2) = (0, 0);
  
  my $read = 0;

  my $samtools = `which samtools`;
  chomp( $samtools);
   
  my $command .= "$samtools view -F4 $infile |  ";
  open ( my $pipe, "$command " ) || die "Could not open '$command': $!\n";
  while(<$pipe>) {
    chomp;
    my (undef, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");    
    
    $read += length($sequence);
    last if ( $sample_size > 0 && $read > $sample_size*2 );

    if ($flags & 0x0080 ) {
      analyse( $sequence, $quality, \%res2);
    }
    else {
      analyse( $sequence, $quality, \%res1);
    }
  }
  close( $pipe );
  
  return( \%res1, \%res2);
}



# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub analyse {
  my ( $seq, $qual, $res) = @_;

  my @seq  = split("", $seq)  if ( $seq );
  my @qual = split("", $qual) if ( $qual );

  my ( $GC, $AT);
  use List::Util qw[max];

  my $length = max (int (@seq || 0), int (@qual || 0));
  
  for(my $i = 0; $i < $length; $i++) {
    if ( $seq ) {
      my $base = uc($seq[$i]);
      $GC++ if ( $base eq 'G' || $base eq 'C');
      $AT++ if ( $base eq 'A' || $base eq 'T');
      $$res{base_dist}[$i]{$base}++;
      $$res{reads}++;
    }

    if ( $qual ) {
      my $base_qual = ord($qual[$i]) - 33;
      $base_qual = int($base_qual/2) if ( $base_qual > 40);
      
      push @{$$res{base_qual}[$i]}, $base_qual;
      $$res{base_qual_dist}{$base_qual}++;
    }
  }

  if ( $AT && $GC  ) {
    my $perc_GC = int( $GC*100/($GC+$AT));
    $perc_GC = int($perc_GC/5) * 5;
    $$res{GC}{ $perc_GC }++ if ( $AT+$GC > 0 );
  }

  return $res;
}



# 
# 
# 
# Kim Brugger (26 Jan 2011), contact: kim.brugger@easih.ac.uk
sub make_plots {
  my ( $data, $title, $outfile ) = @_;
  
  $title ||= "unnamed";
  $outfile ||= "$title";


   if ( $$data{base_dist} ) {
     _plot_base_dist($$data{base_dist}, $$data{reads}, "$outfile.BaseDist");
   }

   if ( $$data{base_qual} ) {
     _plot_base_qual($$data{base_qual}, "$outfile.BaseQual");
   }

   if ( $$data{base_qual_dist} ) {
     _plot_base_qual_dist($$data{base_qual_dist}, "$outfile.BaseQualDist");
   }

  if ( $$data{GC} ) {
    _plot_GC($$data{GC}, "$outfile.GC");
  }

}



# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub _plot_base_qual_dist {
  my ($data, $outfile_prefix) = @_;

  open (my $out, "> $outfile_prefix.R ") || die "Could not open '$outfile_prefix.R': $!\n";
  my (@values, @positions);
  foreach my $QV ( sort { $a <=> $b } keys %$data) {
    push @values, "$QV\t$$data{$QV}\n";
  }

  print $out "\t" . join("\t", @values) . "\n";
  close( $out );  

  open (my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R  "qualHist = read.table('$outfile_prefix.R', header=TRUE, check.names=FALSE)\n";
  print $R "pdf('$outfile_prefix.QualHist.pdf')\n";
  print $R "plot(qualHist, main='$outfile_prefix:  Quality Score Distribution', xlab='Quality score', ylab='Observations', type='h')\n";
  print $R "dev.off()\n";
  close ($R);
  
  system "rm $outfile_prefix.R";
  
}


# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub _plot_base_qual {
  my ($data, $outfile_prefix) = @_;

  open ( my $out, "> $outfile_prefix.R ") || die "Could not open '$outfile_prefix.R': $!\n";
  for(my $i = 0; $i < @{$data}; $i++ ) {
    print $out join("\t", $i+1, @{$$data[$i]}) . "\n";
  }
  close $out;

  open (my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "pdf('$outfile_prefix.pdf')\n";

  print $R "datalines <- readLines('$outfile_prefix.R', n=-1)\n";
  print $R "datalist <- strsplit(datalines, \"\t\")\n";

  print $R "plotlist = list()\n";
  print $R "for (i in 1:length(datalist) ) {\n";
  print $R "plotlist[[datalist[[i]][1]]] <- as.numeric(datalist[[i]][-1])\n";
  print $R "}\n";

  print $R "boxplot(plotlist, ylim=c(0,45), outline=FALSE, main='$outfile_prefix', xlab='QV distribution')\n";
  print $R "dev.off()\n";
  close ($R);

  system "rm $outfile_prefix.R";
}




# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub _plot_base_dist {
  my ($data, $reads, $outfile_prefix) = @_;
  
  open ( my $out, "> $outfile_prefix.R ") || die "Could not open '$outfile_prefix.R': $!\n";

  my (@As, @Cs, @Gs, @Ts, @Ns, @pos);
  
  for(my $i = 0; $i < @{$data}; $i++ ) {
  
    push @pos, $i + 1;
    push @As, ($$data[ $i ]{A} || $$data[ $i ]{0} || 0)/$reads*100;
    push @Cs, ($$data[ $i ]{C} || $$data[ $i ]{1} || 0)/$reads*100;
    push @Gs, ($$data[ $i ]{G} || $$data[ $i ]{2} || 0)/$reads*100;
    push @Ts, ($$data[ $i ]{T} || $$data[ $i ]{3} || 0)/$reads*100;
    push @Ns, ($$data[ $i ]{N} || $$data[ $i ]{'.'} || 0)/$reads*100;
    
  }

  print $out "\t" . join("\t", @pos) . "\n";
  print $out "A\t" . join("\t", @As) . "\n";
  print $out "C\t" . join("\t", @Cs) . "\n";
  print $out "G\t" . join("\t", @Gs) . "\n";
  print $out "T\t" . join("\t", @Ts) . "\n";
  print $out "N\t" . join("\t", @Ns) . "\n";
  
  close ($out);
  
  open (my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "pdf('$outfile_prefix.pdf')\n";
  print $R "baseDist = read.table('$outfile_prefix.R', check.names=FALSE)\n";
  print $R "barplot(as.matrix(baseDist), col=c('green', 'blue', 'black','red', 'grey'), axis.lty=1, main='$outfile_prefix: Base distribution', xlab='Cycle', ylab='Distribution')\n";
  print $R "dev.off()\n";
  close ($R);

  system "rm $outfile_prefix.R";
}


# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub _plot_GC {
  my ($data, $outfile_prefix) = @_;

  my $max_gc = 0;

  open (my $out, "> $outfile_prefix.R ") || die "Could not open '$outfile_prefix.R': $!\n";

  my @values;
  for(my $i = 0; $i <=100; $i+=5) {

    my $value = $$data{ $i } || 0; 
    $max_gc = $value if ( $value > $max_gc );
    
    if ($i - 2.5 < 0 ) {
      push @values, "0\t$value";
    }
    else {
      push @values, ($i - 2.5)."\t$value";
    }

  }
  my $value = $$data{ 100 } || 0; 
  push @values, "100\t$value";
  print $out join("\n", @values) . "\n";

   close( $out );

   $max_gc += 5;

   open ( my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
   print $R "pdf('$outfile_prefix.pdf')\n";
   print $R "GC = read.table('$outfile_prefix.R')\n";
   print $R "plot(spline(GC[c(1,2)]),      type='b',  main='$outfile_prefix: %GC', xlab='%GC distribution', ylab='fraction', col='black', ylim=c(0,$max_gc))\n";
   print $R "dev.off()\n";
   close ($R);

  system "rm $outfile_prefix.R";
}


1;



__END__





