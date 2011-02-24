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
# Does simple QC on a csfasta file without mapping it... If doing a PE run, call the function twice...
# 
# Kim Brugger (03 Aug 2010)
sub csfastaQC  {
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
    $res = analyse( $sequence, undef, $res);
  }


  return $res;
}

# Kim Brugger (03 Aug 2010)
sub qualQC  {
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
    next if (/#/ || /\>/);
    chomp; 

    s/^(\d+)\s*//;
    s/-1/0/g;
    my $qual;
    map { $qual .= chr($_+33)} split(/\s+/);
    
    $read += length( $qual );
    last if ( $sample_size > 0 && $read > $sample_size );

    analyse( undef, $qual, \%res);
    mappable( $qual, \%res); 
  }


  $res{mappability} = sprintf("%.2f", 100*$res{mappable}/$res{quals});

  return \%res;
}
  

# 
# Does simple QC on a fastq file without mapping it... If doing a PE run, call the function twice...
# 
# Kim Brugger (03 Aug 2010)
sub fastQC {
  my ( $infile1, $mappable) = @_;

  $mappable ||= 0;
  my %res;
  my %duplicates;
  my $dup_count = 1;
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
    $duplicates{$sequence}++;
    if ( $read > 1048576*$dup_count) {
      foreach my $seq ( keys %duplicates ) {
	delete $duplicates{$seq} if ($duplicates{$seq} == 1);
      }
      $dup_count++;
    }
    last if ( $sample_size > 0 && $read > $sample_size );
#    print "$read > $sample_size\n";
#    print "$sequence, $strand, $quality\n";
    analyse( $sequence, $quality, \%res);
    mappable( $quality, \%res) if ($mappable); 
  }

  $res{mappability} = sprintf("%.2f", 100*$res{mappable}/$res{reads}) if ($mappable); 
  foreach my $seq ( keys %duplicates ) {
    delete $duplicates{$seq} if ($duplicates{$seq} == 1);
  }
  $res{ duplicates } = \%duplicates;

  return \%res;
}


# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub bamQC {
  my ( $infile, $mappable) = @_;

  my (%res1, %res2);
  my ($read1, $read2) = (0, 0);
  
  $mappable ||= 1;

  my $samtools = `which samtools`;
  chomp( $samtools);
   
  my $command .= "$samtools view -F4 $infile |  ";
  open ( my $pipe, "$command " ) || die "Could not open '$command': $!\n";
  while(<$pipe>) {
    chomp;
    my (undef, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");    
    
    last if ( $sample_size > 0 && $read1 > $sample_size && $read2 > $sample_size );

    if ($flags & 0x0080 ) {
      $read2 += length($sequence);
      analyse( $sequence, $quality, \%res2);
      mappable( $quality, \%res2) if ($mappable); 
    }
    else {
      $read1 += length($sequence);
      analyse( $sequence, $quality, \%res1);
      mappable( $quality, \%res1) if ($mappable); 
    }
  }
  close( $pipe );
  

  $res1{mappability} = sprintf("%.2f", 100*$res1{mappable}/$res1{quals}) if ($mappable);
  $res2{mappability} = sprintf("%.2f", 100*$res2{mappable}/$res2{quals}) if ($mappable);
  return( \%res1, \%res2);
}



# 
# 
# 
# Kim Brugger (28 Jan 2011)
sub mappable {
  my ( $quality, $res) = @_; 
  
  my $cutoff  =  9;
  my $maximum =  11;
#  $maximum = 0.22*length( $quality ) if (length( $quality ) != 50);

  my $count = 0;
  foreach my $QC ( map { (ord($_) - 33) } split("", $quality)) {
    $count++ if ( $QC < $cutoff);
  }
  
  $$res{mappable}++ if ( $count < $maximum);
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
    }

    if ( $qual ) {
      my $base_qual = ord($qual[$i]) - 33;
      $base_qual = int($base_qual/2) if ( $base_qual > 40);
      
      $$res{base_qual}[$i]{$base_qual}++;
      $$res{base_qual_dist}{$base_qual}++;
    }
  }
  $$res{reads}++;

  $$res{quals}++ if ( $qual);

  if ( $AT && $GC  ) {
    my $perc_GC = int( $GC*100/($GC+$AT));
    $perc_GC = int($perc_GC/5) * 5;
    $$res{GC}{ $perc_GC }++ if ( $AT+$GC > 0 );
  }

  if ( $seq ) {
    my $partial_adaptor_mapping = check_for_partial_adaptors($seq);
    my @pam = split("", $partial_adaptor_mapping) ;
    for (my $i=0;$i<@pam;$i++) {
      $$res{partial_adaptor_mapping}{$i} += $pam[$i];
    }
  }

  return $res;
}



# 
# 
# 
# Kim Brugger (26 Jan 2011), contact: kim.brugger@easih.ac.uk
sub make_plots {
  my ( $data, $outfile, $title ) = @_;
  
  $title ||= $outfile;
  $title =~ s/.*\///;

  if ( $$data{base_dist} ) {
    my $ACsplit = _plot_base_dist($$data{base_dist}, $$data{reads}, "$outfile", $title);
    $$data{ACsplit} = $ACsplit;
  }
  
  if ( $$data{base_qual} ) {
    _plot_base_qual($$data{base_qual}, "$outfile", $title);
  }
  
  if ( $$data{base_qual_dist} ) {
    my $Q30 = _plot_base_qual_dist($$data{base_qual_dist}, "$outfile", $title);
    $$data{Q30} = $Q30;
  }
  
  if ( $$data{GC} ) {
    _plot_GC($$data{GC}, "$outfile", $title);
  }
  
  if ( $$data{duplicates} ) {
    my $duplicates = _plot_duplicates($$data{duplicates},  "$outfile", $title);
    $$data{perc_dup} = sprintf("%.2f", 100*$duplicates/$$data{reads});
    foreach my $seq ( keys %{$$data{duplicates}}) {
      my $count = $$data{duplicates}{$seq};
      delete $$data{duplicates}{$seq};
      $$data{duplicates}{$seq}{count}   = $count;
      $$data{duplicates}{$seq}{percent} = sprintf("%.2f", 100*$$data{duplicates}{$seq}{count}/$$data{reads});
    }
  }

  if ($$data{partial_adaptor_mapping}) {
    my $partial_adaptor = _plot_partial_adaptor_mapping($$data{partial_adaptor_mapping}, $$data{reads}, "$outfile", $title);
    $$data{partial_adaptor} = sprintf("%.2f", 100*$partial_adaptor/$$data{reads}) if ($partial_adaptor);
  }
  
}



# 
# 
# 
# Kim Brugger (27 Jan 2011)
sub _plot_partial_adaptor_mapping {
  my ($data, $reads, $outfile_prefix, $title) = @_;

  my $pam_count = 0;

  open (my $out, "> $outfile_prefix\_PAM.R ") || die "Could not open '$outfile_prefix\_PAM.R': $!\n";
  foreach my $pos ( sort { $a <=> $b } keys %$data ) {
    print $out "$pos\t".($$data{$pos}*100/$reads)."\n";
    $pam_count += $$data{$pos};
  }
  close( $out );  

  return if ($pam_count == 0 );

  $pam_count /= int(keys %$data );


  open (my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "PAM = read.table('$outfile_prefix\_PAM.R', header=TRUE, check.names=FALSE)\n";
  print $R "pdf('$outfile_prefix\_PAM.pdf')\n";
  print $R "plot(PAM, main='$title:  Adaptor mapping', xlab='Cycles', ylab='Percent of reads', type='h', ylim=c(0,100))\n";
  print $R "dev.off()\n";
  close ($R);
  
  
  system "rm $outfile_prefix\_PAM.R";

  
  
  return $pam_count;
}


# 
# 
# 
# Kim Brugger (27 Jan 2011)
sub _plot_duplicates {
  my ($data, $outfile_prefix, $title) = @_;

  my %counts;
  my $duplicates = 0;
  foreach my $seq ( keys %$data ) {
    if ($seq =~ /^N*\z/ || $$data{$seq} == 1 ) {
      delete $$data{$seq};
      next;
    }

    $counts{$$data{$seq}}++;
    $duplicates += $$data{$seq};
  }
  
  return $duplicates if ( ! $duplicates );

  open (my $out, "> $outfile_prefix\_DupHist.R ") || die "Could not open '$outfile_prefix\_DupHist.R': $!\n";
  foreach my $count ( sort { $a <=> $b } keys %counts ) {
    print $out "$count\t$counts{$count}\n";
  }
  close( $out );  

  open (my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "dupHist = read.table('$outfile_prefix\_DupHist.R', header=TRUE, check.names=FALSE)\n";
  print $R "pdf('$outfile_prefix\_DupHist.pdf')\n";
  print $R "plot(dupHist, main='$title:  Duplicates Distribution', xlab='Nr of duplications', ylab='Observations', type='h')\n";
  print $R "dev.off()\n";
  close ($R);
  
  system "rm $outfile_prefix\_DupHist.R";
  return $duplicates;
}



# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub _plot_base_qual_dist {
  my ($data, $outfile_prefix, $title) = @_;

  open (my $out, "> $outfile_prefix\_QualHist.R ") || die "Could not open '$outfile_prefix\_QualHist.R': $!\n";
  my (@values, @positions);
  my ( $q30, $total ) = (0,0);
  foreach my $QV ( sort { $a <=> $b } keys %$data) {
    push @values, "$QV\t$$data{$QV}\n";
    $total += $$data{$QV};
    $q30 += $$data{$QV} if ( $QV >= 30 );
  }

  print $out "\t" . join("\t", @values) . "\n";
  close( $out );  

  open (my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R  "qualHist = read.table('$outfile_prefix\_QualHist.R', header=TRUE, check.names=FALSE)\n";
  print $R "pdf('$outfile_prefix\_QualHist.pdf')\n";
  print $R "plot(qualHist, main='$title:  Quality Score Distribution', xlab='Quality score', ylab='Observations', type='h')\n";
  print $R "dev.off()\n";
  close ($R);
  
  system "rm $outfile_prefix\_QualHist.R";

  return sprintf("%.2f", 100*$q30/$total);
}


# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub _plot_base_qual {
  my ($data, $outfile_prefix, $title) = @_;

  open ( my $out, "> $outfile_prefix\_BaseQual.R ") || die "Could not open '$outfile_prefix\_BaseQual.R': $!\n";
  for(my $i = 0; $i < @{$data}; $i++ ) {
    print $out "$i\t";
    foreach my $score ( keys %{$$data[$i]}) {
      for( my $k = 0; $k < $$data[$i]{$score}; $k++) {
	print $out "$score\t";
      }
    }
    print $out "\n";
  }
  close $out;

  open (my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "pdf('$outfile_prefix\_BaseQual.pdf')\n";

  print $R "datalines <- readLines('$outfile_prefix\_BaseQual.R', n=-1)\n";
  print $R "datalist <- strsplit(datalines, \"\t\")\n";

  print $R "plotlist = list()\n";
  print $R "for (i in 1:length(datalist) ) {\n";
  print $R "plotlist[[datalist[[i]][1]]] <- as.numeric(datalist[[i]][-1])\n";
  print $R "}\n";

  print $R "boxplot(plotlist, ylim=c(0,45), outline=FALSE, main='$title: Base qualities', xlab='QV distribution')\n";
  print $R "dev.off()\n";
  close ($R);

  system "rm $outfile_prefix\_BaseQual.R";
}




# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub _plot_base_dist {
  my ($data, $reads, $outfile_prefix, $title) = @_;
  
  open ( my $out, "> $outfile_prefix\_BaseDist.R ") || die "Could not open '$outfile_prefix\_BaseDist.R': $!\n";

  my (@As, @Cs, @Gs, @Ts, @Ns, @pos);

  my $ACsplit = 0;
  my $bases;
  for(my $i = 0; $i < @{$data}; $i++ ) {
  
    push @pos, $i + 1;
    push @As, ($$data[ $i ]{A} || $$data[ $i ]{0} || 0)/$reads*100;
    push @Cs, ($$data[ $i ]{C} || $$data[ $i ]{1} || 0)/$reads*100;
    push @Gs, ($$data[ $i ]{G} || $$data[ $i ]{2} || 0)/$reads*100;
    push @Ts, ($$data[ $i ]{T} || $$data[ $i ]{3} || 0)/$reads*100;
    push @Ns, ($$data[ $i ]{N} || $$data[ $i ]{'.'} || 0)/$reads*100;

    $ACsplit += ($$data[ $i ]{A} || $$data[ $i ]{0} || 0)/$reads*100;
    $ACsplit += ($$data[ $i ]{T} || $$data[ $i ]{3} || 0)/$reads*100;
    $bases++;
  }

  print $out "\t" . join("\t", @pos) . "\n";
  print $out "A\t" . join("\t", @As) . "\n";
  print $out "C\t" . join("\t", @Cs) . "\n";
  print $out "G\t" . join("\t", @Gs) . "\n";
  print $out "T\t" . join("\t", @Ts) . "\n";
  print $out "N\t" . join("\t", @Ns) . "\n";
  
  close ($out);
  
  open (my $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "pdf('$outfile_prefix\_BaseDist.pdf')\n";
  print $R "baseDist = read.table('$outfile_prefix\_BaseDist.R', check.names=FALSE)\n";
  print $R "barplot(as.matrix(baseDist), col=c('green', 'blue', 'black','red', 'grey'), axis.lty=1, main='$title: Base distribution', xlab='Cycle', ylab='Distribution')\n";
  print $R "dev.off()\n";
  close ($R);

  system "rm $outfile_prefix\_BaseDist.R";


  $ACsplit = sprintf("%.2f",$ACsplit/$bases);
  return $ACsplit;
}


# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub _plot_GC {
  my ($data, $outfile_prefix, $title) = @_;

  my $max_gc = 0;

  open (my $out, "> $outfile_prefix\_GC.R ") || die "Could not open '$outfile_prefix\_GC.R': $!\n";

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
   print $R "pdf('$outfile_prefix\_GC.pdf')\n";
   print $R "GC = read.table('$outfile_prefix\_GC.R')\n";
   print $R "plot(spline(GC[c(1,2)]),      type='b',  main='$title: %GC', xlab='%GC distribution', ylab='fraction', col='black', ylim=c(0,$max_gc))\n";
   print $R "dev.off()\n";
   close ($R);

  system "rm $outfile_prefix\_GC.R";
}



my %adaptors = ( 'ACACTCTTTCCCTACACGACGCTGTTCCATCT'                              => 'Illumina Single End Apapter 1',	     
		 'CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT'                            => 'Illumina Single End Apapter 2',	     
		 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'    => 'Illumina Single End PCR Primer 1',     
		 'CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT'                            => 'Illumina Single End PCR Primer 2',     
		 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'  	                   => 'Illumina Single End Sequencing Primer',
		 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'                             => 'Illumina Paired End Adapter 1',  	
		 'CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'                             => 'Illumina Paired End Adapter 2',  	
		 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'    => 'Illumina Paried End PCR Primer 1',  
		 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT' => 'Illumina Paired End PCR Primer 2',  
		 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'                             => 'Illumina Paried End Sequencing Primer 1',  
		 'CGGTCTCGGCATTCCTACTGAACCGCTCTTCCGATCT'                         => 'Illumina Paired End Sequencing Primer 2',  
		 
		 'ACAGGTTCAGAGTTCTACAGTCCGAC'                                    => 'Illumina DpnII expression Adapter 1',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina DpnII expression Adapter 2',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina DpnII expression PCR Primer 1',  
		 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'                  => 'Illumina DpnII expression PCR Primer 2',  
		 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC'                              => 'Illumina DpnII expression Sequencing Primer',		
		 
		 'ACAGGTTCAGAGTTCTACAGTCCGACATG'                                 => 'Illumina NlaIII expression Adapter 1',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina NlaIII expression Adapter 2',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina NlaIII expression PCR Primer 1',  
		 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'                  => 'Illumina NlaIII expression PCR Primer 2',  
		 'CCGACAGGTTCAGAGTTCTACAGTCCGACATG'                              => 'Illumina NlaIII expression Sequencing Primer', 
		 
		 'GTTCAGAGTTCTACAGTCCGACGATC'                                    => 'Illumina Small RNA Adapter 1',  	
		 'TCGTATGCCGTCTTCTGCTTGT'                                        => 'Illumina Small RNA Adapter 2',  	
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina Small RNA RT Primer',  	
		 'CAAGCAGAAGACGGCATACGA'                                         =>  'Illumina Small RNA PCR Primer 1',  
		 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'                  => 'Illumina Small RNA PCR Primer 2',  
		 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC'                              => 'Illumina Small RNA Sequencing Primer',  
		 
		 'GATCGGAAGAGCACACGTCT'                                          => 'Illumina Multiplexing Adapter 1',  
		 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'                             => 'Illumina Multiplexing Adapter 2',  
		 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'    => 'Illumina Multiplexing PCR Primer 1.01',		
		 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'                            => 'Illumina Multiplexing PCR Primer 2.01',		
		 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'                             => 'Illumina Multiplexing Read1 Sequencing Primer',	
		 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'                             => 'Illumina Multiplexing Index Sequencing Primer',	
		 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'                            => 'Illumina Multiplexing Read2 Sequencing Primer',	
		 
		 'CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 1',  	
		 'CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 2',  	
		 'CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 3',  	
		 'CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 4',  	
		 'CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 5',  	
		 'CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 6',  	
		 'CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 7',  	
		 'CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 8',  	
		 'CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 9',  	
		 'CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 10',  	
		 'CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 11',  	
		 'CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC'                   => 'Illumina PCR Primer Index 12',  	
		 
		 'GATCGTCGGACTGTAGAACTCTGAAC'                                    => 'Illumina DpnII Gex Adapter 1',  	
		 'ACAGGTTCAGAGTTCTACAGTCCGAC'                                    => 'Illumina DpnII Gex Adapter 1.01',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina DpnII Gex Adapter 2',  	
		 'TCGTATGCCGTCTTCTGCTTG'                                         => 'Illumina DpnII Gex Adapter 2.01',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina DpnII Gex PCR Primer 1',  
		 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'                  => 'Illumina DpnII Gex PCR Primer 2',  
		 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC'                              => 'Illumina DpnII Gex Sequencing Primer',  
		 
		 'TCGGACTGTAGAACTCTGAAC'                                         => 'Illumina NlaIII Gex Adapter 1.01',  
		 'ACAGGTTCAGAGTTCTACAGTCCGACATG'                                 => 'Illumina NlaIII Gex Adapter 1.02',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina NlaIII Gex Adapter 2.01',  
		 'TCGTATGCCGTCTTCTGCTTG'                                         => 'Illumina NlaIII Gex Adapter 2.02',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina NlaIII Gex PCR Primer 1',  
		 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'                  => 'Illumina NlaIII Gex PCR Primer 2',  
		 'CCGACAGGTTCAGAGTTCTACAGTCCGACATG'                              => 'Illumina NlaIII Gex Sequencing Primer',		
		 
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina Small RNA RT Primer',  	
		 'GTTCAGAGTTCTACAGTCCGACGATC'                                    => 'Illumina 5p RNA Adapter',  	
		 'TCGTATGCCGTCTTCTGCTTGT'                                        => 'Illumina RNA Adapter1',  	
		 
		 'ATCTCGTATGCCGTCTTCTGCTTG'                                      => 'Illumina Small RNA 3p Adapter 1',  
		 'CAAGCAGAAGACGGCATACGA'                                         => 'Illumina Small RNA PCR Primer 1',  
		 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'                  => 'Illumina Small RNA PCR Primer 2',  
		 'CGACAGGTTCAGAGTTCTACAGTCCGACGATC'                              => 'Illumina Small RNA Sequencing Primer',  );


# 
# 
# 
# Kim Brugger (03 Feb 2011)
sub check_for_adaptors {
  my( $seq ) = @_;

  foreach my $adaptor ( %adaptors ) {
    my $rev_adaptor = reverse $adaptor;
    $rev_adaptor =~ tr/[ACGT]/[TGCA]/;
    return $adaptors{ $adaptor } if ( $seq =~ /$adaptor/  || $seq =~ /$rev_adaptor/);
  }

  return "Unkown";
}


# 
# 
# 
# Kim Brugger (03 Feb 2011)
sub check_for_partial_adaptors {
  my( $seq ) = @_;
  
  my $steps       = 5;
  my $frag_length = 15;
  
  my @fragments;
  for(my $i=0;$i<length($seq);$i+=5) {
    push @fragments, substr($seq, $i, $frag_length);
  }


  my $hit_count = 0;
  my $hit_map   = "0" x length($seq);
  for ( my $i = 0; $i<@fragments; $i++ ) {
    my $fragment = $fragments[$i];
    next if (length($fragment ) < $frag_length -2);
    foreach my $adaptor ( %adaptors ) {
      my $rev_adaptor = reverse $adaptor;
      $rev_adaptor =~ tr/[ACGT]/[TGCA]/;
      if ($adaptor =~ /$fragment/ || $rev_adaptor =~ /$fragment/) {
	substr($hit_map, $i*$steps, $frag_length) = '1' x $frag_length;
	$hit_count++;
	goto NXT_FRAG;
	last;
      }
    }
  NXT_FRAG:
  }

#  print "Adaptor hit count: $hit_map, $hit_count\n" if ( $hit_map );
  return $hit_map;
}




1;



__END__





