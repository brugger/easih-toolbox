#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (28 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

use strict;
use Getopt::Long;

my %opts;
getopts('1:2:b:o:hmusp:R:r:f:', \%opts);

my $fq_file       = $opts{f};
my $bam_file      = $opts{b};
my $out_file      = $opts{o};
my $report        = $opts{r} || "";
my $reference     = $opts{R} || die "no reference given (-R)\n";
my $platform      = uc($opts{p}) || die "no platform given (-p[ SOLEXA, SOLID])\n"; 
$platform = 'SOLEXA' if ( $platform eq 'ILLUMINA');

my $limit         = 10000;

my $samtools  = '/usr/local/bin/samtools';

my $ALL_READS     = 0;
my $UNMAPPED_READ = 1;
my $MAPPED_READ   = 2;
my $FIRST_READ    = 3;
my $SECOND_READ   = 4;

my @read_cat = ($ALL_READS, $UNMAPPED_READ, $MAPPED_READ, $FIRST_READ, $SECOND_READ);



my ( %base_qual_dist, @base_dist, @base_qual);
#my ( %mapped_base_qual_dist, @mapped_base_dist, @mapped_base_qual);
#my ( %unmapped_base_qual_dist, @unmapped_base_dist, @unmapped_base_qual);
my ( %GC, %bwa );


my $Q_min = 20;

validate_reference();

#random sample $limit MB out of the fastq file

use File::Temp;
system "mkdir tmp" if ( ! -d './tmp');
my ($tmp_fh, $tmp_file) = File::Temp::tempfile(DIR => "./tmp" );

sample( $fq_file, $tmp_file, 20);
map_and_qc( $tmp_file );

map_and_qc("tmp/PSCn6dr0TF");




#read_bam($bam_file) if ( $bam_file );

die "no input, read the code\n" if ( ! $fq_file && ! $bam_file );



report();

#system "rm -rf tmp";

# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub report {

  my $infile = $fq_file;
  $infile =~ s/.*\///;

  $infile = $report . $infile;

  $out_file ||= $infile;

  my ($Q_min_count) = (0,0);

  my ($out, $R);
  
  my $max_qual = 0;

  open ( $out, "> $out_file.MeanQual ") || die "Could not open '$out_file.MeanQual': $!\n";
  for(my $i = 0; $i < @{$base_qual[$ALL_READS]}; $i++ ) {

    my @values = ($i+1);
    foreach my $category ( @read_cat ) {
      
      if ( ! $bwa{ $category }) {
	push @values, 0;
	next;
      }
      
      if ( $base_qual[$category][ $i ] ) {
	push @values, $base_qual[$category][ $i ] / $bwa{$category};
	$max_qual = $base_qual[$category][ $i ] / $bwa{$category} if ( $max_qual < $base_qual[$category][ $i ] / $bwa{$category});
      }
      else {
	push @values, 0;
      }
    }

    print  join("\t", @values) . "\n";
    print $out join("\t", @values) . "\n";
  }
  close( $out );
  
  $max_qual += 5;

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "meanQual = read.table('$out_file.MeanQual', check.names=FALSE)\n";
#  print $R "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
  print $R "png('$out_file.MeanQual.png')\n";
  print $R "plot(meanQual[,2], main='$infile: Quality By Cycle', type='h', xlab='Cycle', ylab='Mean Quality', ylim=c(0,$max_qual))\n";
  print $R "lines(meanQual[,3], type='l', col='red')\n";
  print $R "lines(meanQual[,4], type='l', col='blue')\n";
  print $R "lines(meanQual[,5], type='l', col='green')\n" if ( $bwa{ $SECOND_READ});
  print $R "lines(meanQual[,6], type='l', col='darkgreen')\n" if ( $bwa{ $SECOND_READ});
  print $R "abline(h=15, col='red')\n";
  print $R "dev.off()\n";
  close ($R);


  open ($out, "> $out_file.QualHist ") || die "Could not open '$out_file.QualHist': $!\n";
  my (@values, @positions);
  foreach my $value ( sort { $a <=> $b } keys %{$base_qual_dist{ $ALL_READS}}) {
    push @positions, $value;

    foreach my $category ( @read_cat ) {
      
      if ( ! $bwa{ $category }) {
	push @{$values[ $category ]}, 0;
	next;
      }
      
      if ( $base_qual_dist{$category}{ $value } ) {
	push @{$values[ $category ]}, $base_qual_dist{$category}{ $value };
      }
      else {
	push @{$values[ $category ]}, 0;
      }
      $Q_min_count++ if ( $base_qual_dist{$category}{ $value } < $Q_min );

    }
  }


  print $out "\t" . join("\t", @positions) . "\n";
  foreach my $category ( @read_cat ) {
    print $out "\t" . join("\t", @{$values[$category]}) . "\n";
  }
  close( $out );  

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R  "qualHist = read.table('$out_file.QualHist', header=TRUE, check.names=FALSE)\n";
  print $R "png('$out_file.QualHist.png')\n";
#  print $R "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
  print $R "barplot(as.matrix(qualHist), col=c('black', 'red', 'blue', 'green', 'darkgreen'), axis.lty=1, main='$infile: Quality distribution', xlab='Quality', ylab='Observations', beside=TRUE)\n";
#  print $R "legend(45, 90, c('All', 'Mapped', 'Unmapped', fill=c('grey', 'blue', 'red')))\n";

#  print $R "plot(qualHist, main='$infile:  Quality Score Distribution', type='h', xlab='Quality score', ylab='Observations')\n";
  print $R "dev.off()\n";
  close ($R);
  
  open ( $out, "> $out_file.BaseDist ") || die "Could not open '$out_file.BaseDist': $!\n";

  my (@As, @Cs, @Gs, @Ts, @Ns, @pos);

  for(my $i = 0; $i < @{$base_dist[$MAPPED_READ]}; $i++ ) {
    my $Ns = $bwa{$MAPPED_READ} - ($base_dist[$MAPPED_READ][ $i ]{A} || 0) - ($base_dist[$MAPPED_READ][ $i ]{C} || 0) - ($base_dist[$MAPPED_READ][ $i ]{G} || 0) - ($base_dist[$MAPPED_READ][ $i ]{T} || 0);

    push @pos, $i + 1;
    push @As, ($base_dist[$MAPPED_READ][ $i ]{A} || 0)/$bwa{$MAPPED_READ}*100;
    push @Cs, ($base_dist[$MAPPED_READ][ $i ]{C} || 0)/$bwa{$MAPPED_READ}*100;
    push @Gs, ($base_dist[$MAPPED_READ][ $i ]{G} || 0)/$bwa{$MAPPED_READ}*100;
    push @Ts, ($base_dist[$MAPPED_READ][ $i ]{T} || 0)/$bwa{$MAPPED_READ}*100;
    push @Ns, ($Ns || 0)/$bwa{$MAPPED_READ}*100;

  }

  print $out "\t" . join("\t", @pos) . "\n";
  print $out "A\t" . join("\t", @As) . "\n";
  print $out "C\t" . join("\t", @Cs) . "\n";
  print $out "G\t" . join("\t", @Gs) . "\n";
  print $out "T\t" . join("\t", @Ts) . "\n";
  print $out "N\t" . join("\t", @Ns) . "\n";

  close ($out);

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "png('$out_file.BaseDist.png')\n";
  print $R "baseDist = read.table('$out_file.BaseDist', check.names=FALSE)\n";
#  print $R "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
  print $R "barplot(as.matrix(baseDist), col=c('green', 'blue', 'black','red', 'grey'), axis.lty=1, main='$infile: Base distribution', xlab='Cycle', ylab='Distribution')\n";
#  print $R "legend(45, 90, rev(rownames(baseDist)), fill=rev(c('green', 'blue', 'black','red', 'grey')))\n";
  print $R "dev.off()\n";
  close ($R);

  my $max_gc = 0;

  open ( $out, "> $out_file.GC ") || die "Could not open '$out_file.GC': $!\n";
  foreach my $value ( sort {$a <=> $b} keys %{$GC{ $ALL_READS}} ) {

    my @values;
    if ($value - 2.5 < 0 ) {
      push @values, "0";
    }
    else {
      push @values, ($value - 2.5);
    }

    foreach my $category ( @read_cat ) {
      
      if ( ! $bwa{ $category }) {
	push @values, 0;
	next;
      }
      
      if ( $GC{$category}{ $value } ) {
#	print "$GC{$category}{ $value }*100 / $bwa{ $category}); ==> " . ($GC{$category}{ $value }*100 / $bwa{ $category}) . "\n";
	push @values, $GC{$category}{ $value }*100 / $bwa{ $category};
	$max_gc = $GC{$category}{ $value }*100/$bwa{ $category} if ( $max_gc < $GC{$category}{ $value }*100/$bwa{ $category});
      }
      else {
	push @values, 0;
      }
    }

    print  join("\t", @values) . "\n";
    print $out join("\t", @values) . "\n";
  }
  close( $out );

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "png('$out_file.GC.png')\n";
  print $R "GC = read.table('$out_file.GC')\n";
  print $R "plot(spline(GC[c(1,2)]),      type='b',  xlab='%GC', ylab='fraction', col='black' ylim=c(0,$max_gc))\n";
  print $R "lines(spline(GC[c(1,3)]),      type='b',  xlab='%GC', ylab='fraction', col='red')\n";
  print $R "lines(spline(GC[c(1,4)]),      type='b',  xlab='%GC', ylab='fraction', col='blue')\n";
  print $R "lines(spline(GC[c(1,5)]),      type='b',  xlab='%GC', ylab='fraction', col='green')\n" if ( $bwa{ $SECOND_READ});
  print $R "lines(spline(GC[c(1,6)]),      type='b',  xlab='%GC', ylab='fraction', col='darkgreen')\n" if ( $bwa{ $SECOND_READ});
  print $R "dev.off()\n";
  close ($R);

  open ( $out, "> $out_file.MapQual ") || die "Could not open '$out_file.MapQual': $!\n";
  foreach my $mapq ( sort {$a <=> $b} keys %{$bwa{mapq}}) {
    print $out "$mapq\t$bwa{mapq}{$mapq}\n";
  }
  close( $out );
    
  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "png('$out_file.MapQual.png')\n";
  print $R "mapQual = read.table('$out_file.MapQual')\n";
  print $R "plot(spline(mapQual),        type='h',  xlab='Mapping Quality', ylab='Fraction', col='black')\n";
  print $R "dev.off()\n";
  close ($R);

  
  if ( $bwa{is} ) {
    
    open ( $out, "> $out_file.IS ") || die "Could not open '$out_file.IS': $!\n";
    print $out join("\n", @{$bwa{is}}) ."\n";
    close( $out );
    
    open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
    print $R "png('$out_file.IS.png')\n";
    print $R "IS = read.table('$out_file.IS')\n";
    print $R "hist(IS,  main='Insert size distribution', xlab='size', ylab='Fraction', col='black')\n";
    print $R "dev.off()\n";
    close ($R);
  }

  open ($out, " > $infile.report ") || die "Could not open '$infile.report': $!\n";
  print $out "$bwa{$ALL_READS} reads\n";;
  print $out "$bwa{ $MAPPED_READ } (".(sprintf("%.2f",$bwa{ $MAPPED_READ }*100/$bwa{$ALL_READS}))."%) mapped\n";
  print $out "$bwa{paired} paired\n"       if ($bwa{paired});
  print $out "$bwa{$FIRST_READ} read1\n"  if ($bwa{$SECOND_READ});
  print $out "$bwa{$SECOND_READ} read2\n" if ($bwa{$SECOND_READ});
  print $out "$bwa{properly_paired} properly paired\n" if ($bwa{properly_paired});
  
  close($out);
}



# 
# 
# 
# Kim Brugger (03 Aug 2010)
sub map_and_qc {
  my ( $infile ) = @_;

  my $bwa      = `which bwa`;
  chomp($bwa);
  my $samtools = `which samtools`;
  chomp( $samtools);
  my $command .= "$bwa aln $reference $infile | $bwa samse $reference - $infile | $samtools view -S - |  ";
  print "$command \n";
  
  open ( my $pipe, "$command " ) || die "Could not open '$command': $!\n";
  while(<$pipe>) {
    chomp;
    my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");    

    $bwa{$ALL_READS}++;

    $bwa{paired}++                 if ($flags & 0x0001);
    $bwa{properly_paired}++        if ($flags & 0x0002 );
    $bwa{mate_unmapped}++          if ($flags & 0x0008 );
    $bwa{not_primary_alignment}++  if ($flags & 0x0100 );
    $bwa{dupliacate}++             if ($flags & 0x0400 );

    if ( $flags & 0x0004 ){
      analyse( $sequence, $quality, $UNMAPPED_READ);
      $bwa{$UNMAPPED_READ}++;
      next;
    }

    $bwa{ $MAPPED_READ }++;
    $bwa{mapq}{$mapq}++      if ( $mapq );
    push @{$bwa{is}}, $insert_size if ( $insert_size ) ;
    
    if ($flags & 0x0080 ) {
      $bwa{$SECOND_READ}++;
      analyse( $sequence, $quality, $SECOND_READ);
    }
    else {
      $bwa{$FIRST_READ}++;
      analyse( $sequence, $quality, $FIRST_READ);
    }
  }
  close( $pipe );
  
}





# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub analyse {
  my ( $seq, $qual, $read_type) = @_;

  $read_type ||= 0;

  my @seq  = split("", $seq);
  my @qual = split("", $qual);

  my ( $GC, $AT);

  for(my $i = 0; $i < @seq; $i++) {
    my $base = uc($seq[$i]);
    my $base_qual = ord($qual[$i]) - 33;

    $GC++ if ( $base eq 'G' || $base eq 'C');
    $AT++ if ( $base eq 'A' || $base eq 'T');

    $base_qual[ $read_type   ][$i] += $base_qual;
    $base_qual[ $ALL_READS   ][$i] += $base_qual;
    $base_qual[ $MAPPED_READ ][$i] += $base_qual if ( $read_type == $FIRST_READ || $read_type == $SECOND_READ);

    $base_qual_dist{ $read_type   }{$base_qual}++;
    $base_qual_dist{ $ALL_READS   }{$base_qual}++;
    $base_qual_dist{ $MAPPED_READ }{$base_qual}++ if ( $read_type == $FIRST_READ || $read_type == $SECOND_READ);

    $base_dist[ $read_type ][$i]{$base}++;
    $base_dist[ $ALL_READS ][$i]{$base}++;
    $base_dist[ $MAPPED_READ ][$i]{$base}++ if ( $read_type == $FIRST_READ || $read_type == $SECOND_READ);
  }

  if ( $AT && $GC  ) {
    
    my $perc_GC = int( $GC*100/($GC+$AT));
    
    $perc_GC = int($perc_GC/5) * 5;

    $GC{ $read_type }{ $perc_GC }++ if ( $AT+$GC > 0 );
    $GC{ $ALL_READS }{ $perc_GC }++;
    $GC{ $MAPPED_READ }{ $perc_GC }++ if ( $read_type == $FIRST_READ || $read_type == $SECOND_READ);
  }
}





sub validate_reference {

  my @errors;

  my @bwa_postfixes = ('amb', 'ann', 'bwt', 'pac', 'rbwt', 'rpac', 'rsa', 'sa');

  push @bwa_postfixes, ( 'nt.amb', 'nt.ann', 'nt.pac')  if ( $platform eq "SOLID");


  # Things related to the reference sequence being used.
  
  my ($dir, $basename, $postfix) = $reference =~ /^(.*)\/(.*?)\.(.*)/;
  foreach my $bwa_postfix ( @bwa_postfixes ) {
    push @errors, "$reference.$bwa_postfix does not exists. Did you run bwa index on $reference?"
	if ( ! -e "$reference.$bwa_postfix");
  }


  push @errors, "Platform must be either SOLEXA or SOLID not '$platform'" if ( $platform ne "SOLEXA" && $platform ne 'SOLID');

  # print the messages and die if critical ones.
  die join("\n", @errors) . "\n"   if ( @errors );
}





# 
# 
# 
# Kim Brugger (03 Aug 2010)
sub sample {
  my ($infile, $outfile, $limit) = @_;

  $limit ||= 1;

  my $compressed = 0;
  if ( $infile =~ /gz\z/) {
    $infile =~ s/.gz//;
    print "infile --> $infile\n";
    print STDERR "Un-compressing the fq file\n";
    $compressed++;
    system "gunzip $infile";
  }

  open (my $out, "> $outfile") || die "Could not open '$outfile' for writing: $!\n";

  my %used;

  my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($infile);

  open (my $file, $infile) || die "Could not open '$infile': $!\n";
  my $read = 0;
  my $goal = $limit*1048576;
  if ( $goal >= $size ) {
    system "cp $infile $outfile";
    return;
  }

  while( $goal > $read ) {
    
    my $random_pos = int(rand($size));
#    print "going to $random_pos\n";
    seek( $file, $random_pos, 0);
    while ( <$file> ) {
      last if ( $_ =~ /^\@/);
    }
    
    if ($_ &&  $_ =~ /^\@\w+:/ && !$used{ $_ }) {
      my $name = $_;
      my $seq  = <$file>;
      my $str  = <$file>;
      my $qual = <$file>;
      
      print $out join("", $name, $seq, $str, $qual);
      $read += length("$name$seq$str$qual");
      $used{$file}++;
    }
    
    
  }
  
  
  if ($compressed) {
    print STDERR "Recompressing the fq file\n";
#    system "gunzip $infile";
  }

}


