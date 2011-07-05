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
use Compress::Zlib;
use IO::Uncompress::Gunzip('$GunzipError');
use IO::Seekable;

use EASIH::QC::db;

my $sample_size   = 0;
my $max_counts    = 0;
my $random_sample = 0;
my $fid           = -1;

# 
# If doing random sampling the number of MB to look at. If -1 or 0 the whole dataset is used 
# 
# Kim Brugger (26 Jan 2011)
sub fid {
  my ($file_id) = @_;
  $fid = $file_id;
  return $fid;
}


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
# If doing random sampling the number of MB to look at. If -1 or 0 the whole dataset is used 
# 
# Kim Brugger (26 Jan 2011)
sub sample_reads {
  my ($counts) = @_;
  $max_counts = $counts if (defined  $counts );

  return $max_counts;
}


# 
# 
# Kim Brugger (26 Jan 2011)
sub random_sample {
  my ($random) = @_;

  $random_sample = $random;

  return $random_sample;
}


# 
# Does simple QC on a csfasta file without mapping it... If doing a PE run, call the function twice...
# 
# Kim Brugger (03 Aug 2010)
sub csfastaQC  {
  my ( $infile1, $res ) = @_;

  if ( $random_sample ) {
    if ( $infile1 =~ /gz/) {
      $res = random_sample_csfasta_gz( $infile1, $res);
    }
    else {
      $res = random_sample_csfasta( $infile1, $res);
    }
  }
  else {
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
  }

  
  return $res;
}

# Kim Brugger (03 Aug 2010)
sub qualQC  {
  my ( $infile1, $res) = @_;
  
  if ( $random_sample ) {
    if ( $infile1 =~ /gz/) {
      $res = random_sample_qual_gz( $infile1, $res);
    }
    else {
      $res = random_sample_qual( $infile1, $res);
    }
  }
  else {
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
      
      $res = analyse( undef, $qual, $res);
      $res = mappable( $qual, $res); 
    }
  }

  $$res{mappability} = sprintf("%.2f", 100*$$res{mappable}/$$res{quals});

  return $res;
}
  

# 
# 
# Kim Brugger (03 Aug 2010)
sub fastQC {
  my ( $infile1, $infile2, $res) = @_;

  my $mappable ||= 0;
  my %duplicates;
  my $dup_count = 1;
  my $read = 0;
  my $file1;
  
  my ($res1, $res2);

  if ( $random_sample ) {
    if ( $infile1 =~ /gz/) {
      ($res1, $res2) = random_sample_fastq_gz_files( $infile1, $infile2);
    }
    else {
      ($res1, $res2) = random_sample_fastq_files( $infile1, $infile2);
    }
  }
  else {
    ($res1, $res2) = sample_fastq_files( $infile1, $infile2);
  }

  my ($QC1, $QC2) = analysis( $res1, $res2 );

  return ($QC1, $QC2);
}



# 
# 
# 
# Kim Brugger (04 Jul 2011)
sub sample_fastq_files {
  my ($infile1, $infile2) = @_;
  
  my ($i1, $i2, @res1, @res2);

  if ($infile1 =~ /gz/) {
    open ($i1, "gunzip -c $infile1 | ") || die "Could not open pipeline: $!\n";
  } 
  else {
    open ($i1, "$infile1 ") || die "Could not open '$infile1': $!\n";
  } 

  if ($infile2 && $infile2 =~ /gz/) {
    open ($i2, "gunzip -c $infile2 | ") || die "Could not open pipeline: $!\n";
  } 
  else {
    open ($i2, "$infile2 ") || die "Could not open '$infile2': $!\n";
  } 

  my $read   = 0;
  my $counts = 0;
  while (<$i1>) { 
    my $name1   = $_;
    my $seq1    = <$i1>;
    my $strand1 = <$i1>;
    my $qual1   = <$i1>;
    $read += length( $seq1 );
    push @res1, [$name1, $seq1, $strand1, $qual1];

    last if ( $sample_size > 0 && $read > $sample_size );
    last if ( $counts++ > $max_counts);
    
    if ( $infile2 ) {
      my $name2   = <$i2>;
      my $seq2    = <$i2>;
      my $strand2 = <$i2>;
      my $qual2   = <$i2>;
      push @res2, [$name2, $seq2, $strand2, $qual2];
    }
  }

  return (\@res1, \@res2);
}



# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub bamQC {
  my ( $infile, $mappable) = @_;

  my ($res1, $res2);
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

    if ($flags & 0x0010 ) {
      $sequence = reverse($sequence);
      $sequence =~ tr/[ACGT]/[TGCA]/;
      $quality = reverse( $quality);
    }

    if ($flags & 0x0080 ) {
      $read2 += length($sequence);
      $res2 = analyse( $sequence, $quality, $res2);
      $res2 = mappable( $quality, $res2) if ($mappable); 
    }
    else {
      $read1 += length($sequence);
      $res1 = analyse( $sequence, $quality, $res1);
      $res1 = mappable( $quality, $res1) if ($mappable); 
    }
  }
  close( $pipe );
  

  $$res1{mappability} = sprintf("%.2f", 100*$$res1{mappable}/$$res1{quals}) if ($mappable);
  $$res2{mappability} = sprintf("%.2f", 100*$$res2{mappable}/$$res2{quals}) if ($mappable);
  return( $res1, $res2);
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
  return $res;
}




# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub analysis {
  my ( $data, $data2 ) = @_;

  use EASIH::Misc;

  my ($res1, $res2);
  my (%duplicates1,%duplicates2);

  my $tmp_file1 = EASIH::Misc::tmp_file();
  my $tmp_file2 = EASIH::Misc::tmp_file();
  open( my $t1, "> $tmp_file1") || die "Could not open '$tmp_file1': $!\n";
  foreach my $entry ( @$data ) {
    
    $duplicates1{ $$entry[1] }++;
    $res1 = analyse( $$entry[1], $$entry[3], $res1);
    $$entry[1] = substr($$entry[1], 0, 36)."\n";
    $$entry[3] = substr($$entry[3], 0, 36)."\n";
    print $t1 join("", @$entry);
  }

  foreach my $seq ( keys %duplicates1 ) {
    delete $duplicates1{$seq} if ($duplicates1{$seq} == 1);
  }
  $$res1{ duplicates } = \%duplicates1;

  if ( $data2 ) {
    open( my $t2, "> $tmp_file2") || die "Could not open '$tmp_file2': $!\n";
    foreach my $entry ( @$data2 ) {
      $$entry[1] = substr($$entry[1], 0, 36)."\n";
      $$entry[3] = substr($$entry[3], 0, 36)."\n";
      print $t2 join("", @$entry);
      $res2 = analyse( $$entry[1], $$entry[2], $res2);
      $duplicates2{ $$entry[1] }++;
    }
    
    foreach my $seq ( keys %duplicates2 ) {
      delete $duplicates2{$seq} if ($duplicates2{$seq} == 1);
    }
    $$res2{ duplicates } = \%duplicates2;
  }


  my $mappings = mappings($tmp_file1) if ( ! $data2 );
  $mappings = mappings($tmp_file1, $tmp_file2) if ( $data2 );

  $$res1{mappings} = $mappings if ($mappings);
  $$res2{mappings} = $mappings if ($mappings && $data2);
  

#  system "rm -f $tmp_file1";
#  system "rm -f $tmp_file2";
#  print "rm -f $tmp_file1\n";
#  print "rm -f $tmp_file2\n";

  return ($res1, $res2);
}



# 
# 
# 
# Kim Brugger (04 Jul 2011)
sub mappings {
  my ($file1, $file2) = @_;

  my %res;

  my %libraries = (
    Human          => "/data/refs/human_1kg/bowtie/human_g1k_v37",
    NCBI_ref       => "/data/refs/archive/human_genomic_transcript/ref_contig",
    NCBI_alt       => "/data/refs/archive/human_genomic_transcript/alt_contig_HuRef",
    NCBI_rna       => "/data/refs/archive/human_genomic_transcript/rna",
    Mouse          => "/data/refs/mm9/bowtie/mm9",
    Rat            => "/data/refs/archive/rn4/bowtie/rn4",
    Ecoli          => "/data/refs/archive/Ecoli/U00096_2",
    MyRDB_bact     => "/data/refs/archive/bacteria/myRDP-bacteria",
    GreenGenes     => "/data/refs/archive/greengenes/greengenes_unaligned",
    Yeast          => "/data/refs/archive/Scerevisiae/yeast",
    Pombe          => "/data/refs/archive/pombe/pombe",
    Viruses        => "/data/refs/archive/viruses/gb_virus_seq",
    PhiX           => "/data/refs/archive/PhiX/PhiX",
    Adapters       => "/data/refs/archive/fastqc_contaminants/Contaminants",
    Vectors        => "/data/refs/archive/UniVec/UniVec",
    Pichia         => "/data/refs/archive/pichia/SO/bowtie/pichia",
    Trypanosoma    => "/data/refs/archive/trypanosoma/Tb927_genome_230210",
    Tryp_gamb      => "/data/refs/archive/trypanosoma/Tbgamb_02_v2",
    Leishmania     => "/data/refs/archive/leishmania/leishmania",
    Ribo_ARB_SSU   => "/data/refs/archive/ribo/hs_ssu_r106_embl",
    Ribo_ARB_LSU   => "/data/refs/archive/ribo/hs_lsu_r106_embl",
    Ribo_UNI       => "/data/refs/archive/ribo/mart_uniprot_ribosomal",
    Ribo_ENS       => "/data/refs/archive/ribo/ens_all_ribo",
    Ribo_GO        => "/data/refs/archive/ribo/amigo_ens_anno_mf",
    MART_RNA_GENES => "/data/refs/archive/ribo/mart_rna_genes",
    MART_RRNA_ONLY => "/data/refs/archive/ribo/mart_rrna_only",
    MART_RRNA_PLUS => "/data/refs/archive/ribo/mart_rrna_mtrrna_pseudorrna",);
  
 
  foreach my $library ( keys %libraries) {
    $library = "Human";
    
    
    my $btie;
    if ($file2) {
      print  "bowtie -S -m 1 --maxins 1000 --chunkmbs 512 $libraries{ $library } -1 $file1 -2 $file2 |\n";
      open ($btie, "bowtie -S -m 1 --maxins 1000 --chunkmbs 512 $libraries{ $library } -1 $file1 -2 $file2 |") or die "Can't launch bowtie: $!";
    }
    else {
      open ($btie, "bowtie -S -m 1 --chunkmbs 512 $libraries{ $library } $file1 |") or die "Can't launch bowtie: $!";
    }
    
    my %dup_counts;
    my ($single, $multiple, $paired) = (0,0,0,0);
    while(<$btie>) {
      next if (/^\@/);
      my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");
      
      next if ( $flags & 0x0004);

      
      # the read is mapped in a proper pair
      if ($flags & 0x0002 ) {
	$paired++;
	my ($p1, $p2) = sort ($pos, $mate_pos);
	$dup_counts{"$p1-$p2"}++;
      }
      
      if ($opts =~ /NM:i:2/) {
	$multiple++;
      }
      else {
	$single++;
      }
    }
    close($btie);

    my $dups = 0;
    map {$dups += $dup_counts{$_} if ( $dup_counts{$_}> 2)} keys %dup_counts;
    $res{ $library }{ single   } = $single;
    $res{ $library }{ multiple } = $multiple;
    $res{ $library }{ paired   } = $paired if ($file2);
    $res{ $library }{ dups     } = $dups   if ($file2);
    last;
  }

  print Dumper( \%res );


  return \%res;


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

    if ( $qual && $qual[$i]) {
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
sub base_dist2db {
  my ($QC_data) = @_;

  return if !($$QC_data{'base_dist'});

  my $data  = $$QC_data{'base_dist'};
  my $reads = $$QC_data{'reads'};
  
  my (@As, @Cs, @Gs, @Ts, @Ns, @pos);

  my $ACsplit = 0;
  my $bases;
  my @res;
  for(my $i = 0; $i < @{$data}; $i++ ) {
    push @res, [$i, ($$data[ $i ]{A} || $$data[ $i ]{0} || 0)/$reads*100,
   		    ($$data[ $i ]{C} || $$data[ $i ]{1} || 0)/$reads*100,
		    ($$data[ $i ]{G} || $$data[ $i ]{2} || 0)/$reads*100,
		    ($$data[ $i ]{T} || $$data[ $i ]{3} || 0)/$reads*100,
		    ($$data[ $i ]{N} || $$data[ $i ]{'.'} || 0)/$reads*100];

    $ACsplit += ($$data[ $i ]{A} || $$data[ $i ]{0} || 0)/$reads*100;
    $ACsplit += ($$data[ $i ]{T} || $$data[ $i ]{3} || 0)/$reads*100;
    $bases++;
  }


  EASIH::QC::db::add_basedists($fid, \@res);

  $ACsplit = sprintf("%.2f",$ACsplit/$bases);
  EASIH::QC::db::update_file($fid, undef, $$QC_data{reads}, undef, undef, undef, $ACsplit);
}


# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub base_qual2db {
  my ($QC_data) = @_;

  use Statistics::Descriptive;

  my @res;
  for(my $x = 0; $x < @{$$QC_data{'base_qual'}}; $x++ ) {
    my @counts;
    foreach my $score ( keys %{$$QC_data{'base_qual'}[$x]}) {
      for( my $k = 0; $k < $$QC_data{'base_qual'}[$x]{$score}; $k++) {
	push @counts, $score;
      }
    }

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@counts);
    push @res, [$x, 
		$stat->quantile(0),
		$stat->quantile(1),
		$stat->quantile(2),
		$stat->quantile(3),
		$stat->quantile(4)];
  }

  EASIH::QC::db::add_qvs_boxplot( $fid, \@res );

}


# 
# 
# 
# Kim Brugger (26 Jan 2011)
sub mappings2db {
  my ($QC_data, $fid2) = @_;

  print Dumper($$QC_data{'mappings'});

  foreach my $db (keys %{$$QC_data{mappings}}) {
    EASIH::QC::db::add_mapping_stats( $fid, $fid2, $db, $$QC_data{mappings}{$db}{single}, $$QC_data{mappings}{$db}{multiple}, $$QC_data{mappings}{$db}{dups} );
  }

}




# 
# 
# 
# Kim Brugger (01 Jul 2011)
sub base_qual_dist2db {
  my ($QC_data) = @_;

  my @res;
  my ( $q30, $total ) = (0, 0);
  foreach my $x ( keys %{$$QC_data{'base_qual_dist'}}) {
    push @res, [$x, $$QC_data{base_qual_dist}{ $x }];
    $q30   += $$QC_data{'base_qual_dist'}{$x} if ( $x >= 30 );
    $total += $$QC_data{'base_qual_dist'}{$x};
  }

  EASIH::QC::db::add_qvs_histogram( $fid, \@res );
  $q30 = sprintf("%.2f", 100*$q30/$total);
  EASIH::QC::db::update_file($fid, undef, $$QC_data{reads}, sprintf("%2.f", $q30), undef, undef, undef);
 
}



# 
# 
# 
# Kim Brugger (01 Jul 2011)
sub GC2db {
  my ($QC_data) = @_;

  my $bucket = 5;
  
  my %binned;
  my $total = 0;
  foreach my $x (keys %{$$QC_data{ GC }} ) {
    $binned{ int( $x/$bucket)} += $$QC_data{ GC }{$x};
    $total += $$QC_data{ GC }{$x};
  }

  my @res;
  foreach my $x (keys %binned ) {
    push @res, [$x, sprintf("%.2f", 100*$binned{ $x }/$total)];
  }

  EASIH::QC::db::add_gc_distribution( $fid, \@res );

}



# 
# 
# 
# Kim Brugger (01 Jul 2011)
sub duplicates2db {
  my ($QC_data) = @_;

  my %counts;
  my %dup_seqs;
  my $duplicates = 0;
  foreach my $seq ( keys %{$$QC_data{duplicates}} ) {
    if ($seq =~ /^N*\z/ || $$QC_data{duplicates}{$seq} == 1 ) {
      delete $$QC_data{duplicates}{$seq};
      next;
    }

    $counts{$$QC_data{duplicates}{$seq}}++;
    $dup_seqs{ $seq }++;
    $duplicates += $$QC_data{duplicates}{$seq};
  }

  my @res;
  foreach my $seq (  keys %dup_seqs ) {
    push @res, [ $seq, sprintf("%.2f", 100*$dup_seqs{ $seq }/ $$QC_data{reads}), check_for_adaptors( $seq )];
  }

  EASIH::QC::db::add_duplicate_seqs( $fid, \@res);
  
  my $perc_dup = sprintf("%.2f", 100*$duplicates/($$QC_data{reads}+1));
  EASIH::QC::db::update_file($fid, undef, $$QC_data{reads}, undef, $perc_dup, undef, undef);


  @res = ();
  foreach my $count (  keys %counts ) {
    push @res, [ $count, $counts{$count}];
  }
  
  EASIH::QC::db::add_duplicates( $fid, \@res);
  
}


# 
# 
# 
# Kim Brugger (01 Jul 2011)
sub partial_adaptor2db {
  my ($QC_data) = @_;

  my $pam_count = 0;

  print Dumper($$QC_data{partial_adaptor_mapping});

  my @res;
  foreach my $pos ( keys %{$$QC_data{partial_adaptor_mapping}} ) {
    push @res, [$pos, sprintf("%.2f", $$QC_data{partial_adaptor_mapping}{$pos}*100/$$QC_data{reads})];
    $pam_count += $$QC_data{partial_adaptor_mapping}{$pos};
  }

  print Dumper(\@res);

  EASIH::QC::db::add_adaptors( $fid, \@res);
  $pam_count /= int(keys %{$$QC_data{partial_adaptor_mapping}} );
  $pam_count = sprintf("%.2f", $pam_count*100/$$QC_data{reads});
  print "$pam_count\n";
  EASIH::QC::db::update_file($fid, undef, $$QC_data{reads}, undef, undef, $pam_count, undef);
 
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



#
#
#
# Kim Brugger (03 Aug 2010)
sub random_sample_fastq_gz_files {
  my ($infile1, $infile2 ) = @_;

  my $z1 = new IO::Uncompress::Gunzip($infile1) or die "gunzip failed: $GunzipError\n";
  my $z2 = new IO::Uncompress::Gunzip($infile2) or die "gunzip failed: $GunzipError\n" if ( $infile2 );

  my $read   = 0;
  my $counts = 0;
  my $random_pos = int(rand(1048576));
  
  my (@res1, @res2);


  while( $sample_size > $read && $max_counts > $counts++ ) {
    
    $random_pos = $z1->tell() + int(rand(5120));
    if ( ! $z1->seek($random_pos, SEEK_SET) ) {
      close($z1);
      close($z2) if ( $infile2 );
      $z1 = new IO::Uncompress::Gunzip($infile1) or die "gunzip failed: $GunzipError\n";
      $z2 = new IO::Uncompress::Gunzip($infile2) or die "gunzip failed: $GunzipError\n" if ( $infile2 );
      $random_pos = int(rand(1048576));
      next;
    }
    $z2->seek($random_pos, SEEK_SET) if ( $infile2 );

    while ( $_ =  $z1->getline() ) {
      last if ( $_ =~ /^\@[A-Za-z0-9-_]*:\d+:/ ||  $_ =~ /^\@[A-Za-z0-9-_]*:\d+_\d+_\d+/);
      $z2->getline();
    }

    next if (! $_);

    if ($_ && 
	($_ =~ /^\@[A-Za-z0-9-_]*:\d+:/ ||  
	 $_ =~ /^\@[A-Za-z0-9-_]*:\d+_\d+_\d+/)) {
      my $name1 = $_;
      my $seq1  = $z1->getline();
      my $str1  = $z1->getline();
      my $qual1 = $z1->getline();
      
      push @res1, [$name1, $seq1, $str1, $qual1];
      $read += length($seq1);

      if ( $infile2 ) {
	my $name2 = $z2->getline();;
	my $seq2  = $z2->getline();
	my $str2  = $z2->getline();
	my $qual2 = $z2->getline();
	push @res2, [$name2, $seq2, $str2, $qual2];
      }
    }
  }

  close($z1);
  close($z2) if ( $infile2 );

  return (\@res1, \@res2);
}


sub random_sample_fastq_gz {
  my ($infile, $res ) = @_;

  my $z = new IO::Uncompress::Gunzip($infile) or die "gunzip failed: $GunzipError\n";

  my $read = 0;
  my $random_pos = int(rand(1048576));
  
  while( $sample_size > $read ) {

    
    $random_pos = $z->tell() + int(rand(5120));
    if ( ! $z->seek($random_pos, SEEK_SET) ) {
      close($z);
      $z = new IO::Uncompress::Gunzip($infile) or die "gunzip failed: $GunzipError\n";
      $random_pos = int(rand(1048576));
      next;
    }

    while ( $_ =  $z->getline() ) {
      last if ( $_ =~ /^\@[A-Za-z0-9-_]*:\d+:/ ||  $_ =~ /^\@[A-Za-z0-9-_]*:\d+_\d+_\d+/);
    }

    next if (! $_);

    if ($_ && 
        ($_ =~ /^\@[A-Za-z0-9-_]*:\d+:/ ||  
         $_ =~ /^\@[A-Za-z0-9-_]*:\d+_\d+_\d+/)) {
      my $name = $_;
      my $seq  = $z->getline();
      my $str  = $z->getline();
      my $qual = $z->getline();

      $res = analyse( $seq, $qual, $res);
      $res = mappable( $qual, $res);
      
      $read += length($seq);
    }
  }


  return $res;
}



#
#
#
# Kim Brugger (03 Aug 2010)
sub random_sample_fastq_files {
  my ($infile1, $infile2) = @_;

  my (@res1, @res2);

  my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($infile1);
  
  open (my $file1, $infile1) || die "Could not open '$infile1': $!\n";
  open (my $file2, $infile2) || die "Could not open '$infile2': $!\n" if ( $infile2);
  my $read = 0;
  my $counts = 0;
  while( $sample_size > $read  && $max_counts > $counts) {


    my $random_pos = int(rand($size));
    seek( $file1, $random_pos, 0);
    seek( $file2, $random_pos, 0) if ( $infile2 );
    while ( <$file1> ) {
      last if ( $_ =~ /^\@[A-Za-z0-9-_]*:\d+:/ ||  $_ =~ /^\@[A-Za-z0-9-_]*:\d+_\d+_\d+/);
      <$file2> if ( $infile2 );
    }
    
    if ($_ && ($_ =~ /^\@[A-Za-z0-9-_]*:\d+:/ || $_ =~ /^\@[A-Za-z0-9-_]*:\d+_\d+_\d+/)) {
      my $name1 = $_;
      my $seq1  = <$file1>;
      my $str1  = <$file1>;
      my $qual1 = <$file1>;
      
      push @res1, [$name1, $seq1, $str1, $qual1];

      if ( $infile2 ) {

	my $name2 = <$file2>;
	my $seq2  = <$file2>;
	my $str2  = <$file2>;
	my $qual2 = <$file2>;

	push @res2, [$name2, $seq2, $str2, $qual2];
      }
      
      $read += length($seq1);
    }
  }


  return (\@res1, \@res2);
}



#
#
#
# Kim Brugger (03 Aug 2010)
sub random_sample_csfasta_gz {
  my ($infile, $res ) = @_;

  my $z = new IO::Uncompress::Gunzip($infile) or die "gunzip failed: $GunzipError\n";

  my $read = 0;
  my $random_pos = int(rand(1048576));
  
  while( $sample_size > $read ) {
    
    $random_pos = $z->tell() + int(rand(5120));
    if ( ! $z->seek($random_pos, SEEK_SET) ) {
      close($z);
      $z = new IO::Uncompress::Gunzip($infile) or die "gunzip failed: $GunzipError\n";
#      print STDERR "Start from the beginning again\n";
      $random_pos = int(rand(1048576));
      next;
    }

    while ( $_ =  $z->getline() ) {
      last if ( $_ =~ /^\>/);
    }
    
    if ($_ && $_ =~ /^\>/) {
      my $seq  = $z->getline();
      $seq = substr($seq, 1);

      $res = analyse( $seq, undef, $res);
      
      $read += length($seq);
    }
  }

  return $res;
}


#
#
#
# Kim Brugger (03 Aug 2010)
sub random_sample_csfasta {
  my ($infile1, $res) = @_;

  my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($infile1);
  
  open (my $file1, $infile1) || die "Could not open '$infile1': $!\n";
  my $read = 0;
  while( $sample_size > $read ) {

    my $random_pos = int(rand($size));
    seek( $file1, $random_pos, 0);
    
    while ( <$file1> ) {
      last if ( $_ =~ /^\>/);
    }
    
    if ($_ && $_ =~ /^>/) {
      my $seq  = <$file1>;
      $seq = substr($seq, 1);
      $res = analyse( $seq, undef, $res);
      
      $read += length($seq);
    }
  }


  return $res;
}



#
#
#
# Kim Brugger (03 Aug 2010)
sub random_sample_qual_gz {
  my ($infile, $res ) = @_;

  my $z = new IO::Uncompress::Gunzip($infile) or die "gunzip failed: $GunzipError\n";

  my $read = 0;
  my $random_pos = int(rand(1048576));
  
  while( $sample_size > $read ) {
    
    $random_pos = $z->tell() + int(rand(5120));
    if ( ! $z->seek($random_pos, SEEK_SET) ) {
      close($z);
      $z = new IO::Uncompress::Gunzip($infile) or die "gunzip failed: $GunzipError\n";
      #print STDERR "Start from the beginning again\n"; 
      $random_pos = int(rand(1048576));
      next;
    }

    while ( $_ =  $z->getline() ) {
      last if ( $_ =~ /^\>/);
    }
    
    if ($_ && $_ =~ /^\>/) {
      $_  = $z->getline();
      chomp;
      
      s/^(\d+)\s*//;
      s/-1/0/g;
      my $qual;
      map { $qual .= chr($_+33)} split(/\s+/, $_);
    
      $read += length( $qual );
      last if ( $sample_size > 0 && $read > $sample_size );

      $res = analyse( undef, $qual, $res);
      $res = mappable( $qual, $res); 
      $read += length($qual);
    }
  }

  return $res;
}


#
#
#
# Kim Brugger (03 Aug 2010)
sub random_sample_qual {
  my ($infile1, $res) = @_;

  my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($infile1);
  
  open (my $file1, $infile1) || die "Could not open '$infile1': $!\n";
  my $read = 0;
  while( $sample_size > $read ) {

    my $random_pos = int(rand($size));
    seek( $file1, $random_pos, 0);
    
    while ( <$file1> ) {
      last if ( $_ =~ /^\>/);
    }
    
    if ($_ && $_ =~ /^>/) {
      $_ = <$file1>;
      chomp;
      
      s/^(\d+)\s*//;
      s/-1/0/g;
      my $qual;
      map { $qual .= chr($_+33)} split(/\s+/, $_);
    
      $read += length( $qual );
      last if ( $sample_size > 0 && $read > $sample_size );

      $res = analyse( undef, $qual, $res);
      $res = mappable( $qual, $res); 
      $read += length($qual);
    }
  }

  return $res;
}

1;



__END__







