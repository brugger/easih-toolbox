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

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::QC;

my %opts;
getopts('b:c:f:q:hs:', \%opts);

my $fastq_file    = $opts{f};
my $csfastq_file  = $opts{c};
my $qual_file     = $opts{q};
my $bam_file      = $opts{b};
my $sample_size   = $opts{s} || 10;

#random sample $limit MB out of the fastq file
EASIH::QC::sample_size($sample_size);

if ( $fastq_file ) {
  my $QC = EASIH::QC::fastQC( $fastq_file );
  $fastq_file =~ s/(.*?)\..*/$1/;
  EASIH::QC::make_plots($QC, $fastq_file);
}
if ( $csfastq_file ) {
  my $QC = EASIH::QC::csfastQC( $csfastq_file);
  $csfastq_file =~ s/(.*?)\..*/$1/;
  EASIH::QC::make_plots($QC, $csfastq_file);
}
if ( $qual_file ) {
  my $QC = EASIH::QC::qualQC( $qual_file );
  $qual_file =~ s/(.*?)\..*/$1/;
  EASIH::QC::make_plots($QC, $qual_file);
}
if ( $bam_file ) {
  my ($QC1, $QC2) = EASIH::QC::bamQC( $bam_file );
  $bam_file =~ s/(.*?)\..*/$1/;
  EASIH::QC::make_plots($QC1, "$bam_file.1");
  EASIH::QC::make_plots($QC2, "$bam_file.2");
}


__END__

my $second_file   = $opts{2};
my $bam_file      = $opts{b};
my $no_mapping    = $opts{n} ||  0;
my $sample_size   = $opts{s} || 50; # This is in MB
my $out_file      = $opts{o};
my $report        = $opts{r} || "";
my $reference     = $opts{R} || usage() if ( ! $bam_file && ! $no_mapping);
my $platform      = uc($opts{p}) || die "no platform given (-p[ SOLEXA, SOLID])\n"; 
$platform = 'SOLEXA' if ( $platform eq 'ILLUMINA');

my $samtools  = '/usr/local/bin/samtools';

usage() if ( ! $first_file && ! $bam_file);

validate_input();

use File::Temp;
system "mkdir tmp" if ( ! -d './tmp');
my ($tmp_fh, $tmp_file) = File::Temp::tempfile(DIR => "./tmp" );


if ( $bam_file ) {
  qc_premapped( $bam_file, $sample_size );
}
elsif ( $first_file )  {
  print "Start sampling ..\n";
  sample( $first_file, $second_file, "$tmp_file.1", "$tmp_file.2", $sample_size);
  print "done.....\n";
#  qc( "$first_file", "$second_file" ) if ( $no_mapping );
  qc( "$tmp_file.1", "$tmp_file.2" ) if ( $no_mapping );
  map_and_qc( "$tmp_file.1", "$tmp_file.2" ) if ( !$no_mapping );;
}


report();

#system "rm -rf tmp";

# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub report {

  my $infile = $first_file || $bam_file;
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

    print $out join("\t", @values) . "\n";
  }
  close( $out );
  
  $max_qual += 5;

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "meanQual = read.table('$out_file.MeanQual', check.names=FALSE)\n";
#  print $R "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
  print $R "png('$out_file.MeanQual.png')\n";
  print $R "plot(meanQual[,2], main='$out_file: Quality By Cycle', type='h', xlab='Cycle', ylab='Mean Quality', ylim=c(0,$max_qual))\n";
  print $R "lines(meanQual[,3], type='l', col='red')\n" if ( $bwa{ $UNMAPPED_READ});
  print $R "lines(meanQual[,4], type='l', col='blue')\n"if ( $bwa{ $UNMAPPED_READ} && ! $bwa{ $SECOND_READ});
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

      $Q_min_count++ if ( $base_qual_dist{$category}{ $value } && $base_qual_dist{$category}{ $value } < $Q_min );
    }
  }


  print $out "\t" . join("\t", @positions) . "\n";
  foreach my $category ( @read_cat ) {
    next if ( ! $bwa{ $category });
    print $out "\t" . join("\t", @{$values[$category]}) . "\n";
  }
  close( $out );  

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R  "qualHist = read.table('$out_file.QualHist', header=TRUE, check.names=FALSE)\n";
  print $R "png('$out_file.QualHist.png')\n";
#  print $R "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
  print $R "barplot(as.matrix(qualHist), col=c('black', 'red', 'blue', 'green', 'darkgreen'), axis.lty=1, main='$out_file: Quality distribution', xlab='Quality', ylab='Observations', beside=TRUE)\n" if ( $bwa{ $UNMAPPED_READ} && $bwa{ $SECOND_READ});
  print $R "barplot(as.matrix(qualHist), col=c('black', 'green', 'darkgreen'), axis.lty=1, main='$out_file: Quality distribution', xlab='Quality', ylab='Observations', beside=TRUE)\n" if ( !$bwa{ $UNMAPPED_READ} && $bwa{ $SECOND_READ});
  print $R "barplot(as.matrix(qualHist), col=c('black'), axis.lty=1, main='$out_file: Quality distribution', xlab='Quality', ylab='Observations', beside=TRUE)\n" if ( !$bwa{ $UNMAPPED_READ} && ! $bwa{ $SECOND_READ});
#  print $R "legend(45, 90, c('All', 'Mapped', 'Unmapped', fill=c('grey', 'blue', 'red')))\n";

#  print $R "plot(qualHist, main='$out_file:  Quality Score Distribution', type='h', xlab='Quality score', ylab='Observations')\n";
  print $R "dev.off()\n";
  close ($R);
  
  open ( $out, "> $out_file.BaseDist ") || die "Could not open '$out_file.BaseDist': $!\n";

  my (@As, @Cs, @Gs, @Ts, @Ns, @pos);


  if ( $bwa{$MAPPED_READ} ) {
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
  print $R "barplot(as.matrix(baseDist), col=c('green', 'blue', 'black','red', 'grey'), axis.lty=1, main='$out_file: Base distribution', xlab='Cycle', ylab='Distribution')\n";
#  print $R "legend(45, 90, rev(rownames(baseDist)), fill=rev(c('green', 'blue', 'black','red', 'grey')))\n";
  print $R "dev.off()\n";
  close ($R);
  }

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

    print $out join("\t", @values) . "\n";
  }
  close( $out );

  $max_gc += 5;

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "png('$out_file.GC.png')\n";
  print $R "GC = read.table('$out_file.GC')\n";
  print $R "plot(spline(GC[c(1,2)]),      type='b',  main='$out_file: %GC', xlab='%GC distribution', ylab='fraction', col='black', ylim=c(0,$max_gc))\n";
  print $R "lines(spline(GC[c(1,3)]),      type='b',  xlab='%GC', ylab='fraction', col='red')\n" if ( $bwa{ $UNMAPPED_READ});
  print $R "lines(spline(GC[c(1,4)]),      type='b',  xlab='%GC', ylab='fraction', col='blue')\n" if ( $bwa{ $UNMAPPED_READ} && ! $bwa{ $SECOND_READ});
  print $R "lines(spline(GC[c(1,5)]),      type='b',  xlab='%GC', ylab='fraction', col='green')\n" if ( $bwa{ $SECOND_READ});
  print $R "lines(spline(GC[c(1,6)]),      type='b',  xlab='%GC', ylab='fraction', col='darkgreen')\n" if ( $bwa{ $SECOND_READ});
  print $R "dev.off()\n";
  close ($R);

  if ( $bwa{mapq}) {
    open ( $out, "> $out_file.MapQual ") || die "Could not open '$out_file.MapQual': $!\n";
    foreach my $mapq ( sort {$a <=> $b} keys %{$bwa{mapq}}) {
      print $out "$mapq\t$bwa{mapq}{$mapq}\n";
    }
    close( $out );
    
    open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
    print $R "png('$out_file.MapQual.png')\n";
    print $R "mapQual = read.table('$out_file.MapQual')\n";
    print $R "plot(spline(mapQual),        type='h',  main='$out_file: Mapping quality distribution', xlab='Quality score', ylab='count', col='black')\n";
    print $R "dev.off()\n";
    close ($R);
  }
  
  if ( $bwa{is} ) {
    
    open ( $out, "> $out_file.IS ") || die "Could not open '$out_file.IS': $!\n";
    print $out join("\n", @{$bwa{is}}) ."\n";
    close( $out );
    
    open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
    print $R "png('$out_file.IS.png')\n";
    print $R "IS = read.table('$out_file.IS')\n";
    print $R "hist(IS[,1],  main='Insert size distribution', xlab='size', ylab='Fraction', breaks=300)\n";
    print $R "dev.off()\n";
    close ($R);
  }
  
  print Dumper(\%bwa);

  open ($out, " > $out_file.report ") || die "Could not open '$out_file.report': $!\n";
  print $out "$bwa{$ALL_READS} reads\n";;
  print $out "$bwa{ $MAPPED_READ } (".(sprintf("%.2f",$bwa{ $MAPPED_READ }*100/$bwa{$ALL_READS}))."%) mapped\n" if ($bwa{ $MAPPED_READ });
  print $out "$bwa{paired} paired\n"       if ($bwa{paired});
  print $out "$bwa{$FIRST_READ} read1\n"  if ($bwa{$SECOND_READ});
  print $out "$bwa{$SECOND_READ} read2\n" if ($bwa{$SECOND_READ});
  print $out "$bwa{properly_paired} properly paired\n" if ($bwa{properly_paired});
  print $out "$bwa{duplicate} duplicated reads\n" if ($bwa{duplicate});
  
  close($out);
}





sub validate_input {

  my @errors;

  my @bwa_postfixes = ('amb', 'ann', 'bwt', 'pac', 'rbwt', 'rpac', 'rsa', 'sa');

  push @bwa_postfixes, ( 'nt.amb', 'nt.ann', 'nt.pac')  if ( $platform eq "SOLID");


  # Things related to the reference sequence being used.
  if ( $reference ) {
    my ($dir, $basename, $postfix) = $reference =~ /^(.*)\/(.*?)\.(.*)/;
    foreach my $bwa_postfix ( @bwa_postfixes ) {
      push @errors, "$reference.$bwa_postfix does not exists. Did you run bwa index on $reference?"
	  if ( ! -e "$reference.$bwa_postfix");
    }
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
  my ($infile1, $infile2, $outfile1, $outfile2, $limit) = @_;

  $limit ||= 1;

  if ( $infile1 =~ /gz\z/) {
    print "infile --> $infile1\n";
    my $outfile = $infile1;
    $outfile =~ s/.gz//;
    print STDERR "Un-compressing the fq file $infile1 -> tmp/$outfile\n";
    system "gunzip -c $infile1 > tmp/$outfile";
    $infile1 = "tmp/$outfile";
  }

  if ( $infile2 && $infile2 =~ /gz\z/) {
    print "infile --> $infile2\n";
    my $outfile = $infile2;
    $outfile =~ s/.gz//;
    print STDERR "Un-compressing the fq file $infile2 -> tmp/$outfile\n";
    system "gunzip -c $infile2 > tmp/$outfile";
    $infile2 = "tmp/$outfile";
  }

  open (my $out1, "> $outfile1") || die "Could not open '$outfile1' for writing: $!\n";
  open (my $out2, "> $outfile2") || die "Could not open '$outfile2' for writing: $!\n" if ( $infile2);
  
  my %used;
  
  my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($infile1);
  
  open (my $file1, $infile1) || die "Could not open '$infile1': $!\n";
  open (my $file2, $infile2) || die "Could not open '$infile2': $!\n" if ( $infile2 );
  my $read = 0;
  my $goal = $limit*1048576;
  if ( $goal >= $size ) {
    system "cp $infile1 $outfile1";
    system "cp $infile2 $outfile2" if ( $infile2);
    return;
  }
  
  
  while( $goal > $read ) {
    

    my ($name1, $name2);
    
    my $random_pos = int(rand($size));
#    print "going to $random_pos $goal vs $read\n";

    seek( $file1, $random_pos, 0);
#    print "Finished the jump, now scanning...\n";
    
    while ( <$file1> ) {
#      print ;
      last if ( $_ =~ /^\@[A-Za-z0-9-_]*:\d+:/);
    }
#    print "Found it\n";
    
    if ($_ && $_ =~ /^\@[A-Za-z0-9-_]*:\d+:/ && !$used{ $_ }) {
      my $name = $_;
      my $seq = <$file1>;
      my $str = <$file1>;
      my $qual = <$file1>;
      
      print $out1 join("", $name, $seq, $str, $qual);
      $read += length("$name$seq$str$qual");
      $used{$infile1}++;
      $name1 = $name;
    }

    if ( $infile2 && $name1) {

      die "I should not be here \n";
      
      $name1 =~ s/\/\d\n//;
      
      my %seen;
      
    REREAD:
      
      seek( $file2, $random_pos, 0);
      while ( <$file2> ) {
	last if ( $_ =~ /^\@\w+:\d+:/);
      }
      
      if ($_ && $_ =~ /^\@\w+:\d+:/ && !$used{ $_ }) {
	my $name = $_;
	$name2 = $name;
	$name2 =~ s/\/\d\n//;
	
	next if ($seen{$name2} && $seen{$name2} > 10);
	
	
	if ( $name1 eq $name2 ) {
 print "$name1 eq $name2\n";
	  my $seq = <$file2>;
	  my $str = <$file2>;
	  my $qual = <$file2>;
	  
	  print $out2 join("", $name, $seq, $str, $qual);
	  $read += length("$name$seq$str$qual");
	}
	elsif ($name1 lt $name2 ) {
 print "$name1 lt $name2\n";
	  $random_pos -= 100;
	  $seen{ $name2 }++;
	  goto REREAD;
	}
	else {
 print "$name1 gt $name2\n";
	  $random_pos = tell($file2);
	  $seen{ $name2 }++;
	  goto REREAD;
	}
      }
    }
  }
}


# 
# 
# 
# Kim Brugger (03 Aug 2010)
sub sample_old {
  my ($infile1, $infile2, $outfile1, $outfile2, $limit) = @_;

  $limit ||= 1;

  if ( $infile1 =~ /gz\z/) {
    print "infile --> $infile1\n";
    my $outfile = $infile1;
    $outfile =~ s/.gz//;
    $outfile =~ s/.*\///;
    print STDERR "Un-compressing the fq file $infile1 -> tmp/$outfile\n";
    system "gunzip -c $infile1 > tmp/$outfile";
    $infile1 = "tmp/$outfile";
  }

  if ( $infile2 && $infile2 =~ /gz\z/) {
    print "infile --> $infile2\n";
    my $outfile = $infile2;
    $outfile =~ s/.gz//;
    print STDERR "Un-compressing the fq file $infile2 -> tmp/$outfile\n";
    system "gunzip -c $infile2 > tmp/$outfile";
    $infile2 = "tmp/$outfile";
  }

  open (my $out1, "> $outfile1") || die "Could not open '$outfile1' for writing: $!\n";
  open (my $out2, "> $outfile2") || die "Could not open '$outfile2' for writing: $!\n"   if ( $infile2);

  my %used;

  my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($infile1);

  open (my $file1, $infile1) || die "Could not open '$infile1': $!\n";
  open (my $file2, $infile2) || die "Could not open '$infile2': $!\n" if ( $infile2 );
  my $read = 0;
  my $goal = $limit*1048576;
  if ( $goal >= $size ) {
    system "cp $infile1 $outfile1";
    system "cp $infile2 $outfile2" if ( $infile2);
    return;
  }


  while( $goal > $read ) {

    my ($name1, $name2);
    
    my $random_pos = int(rand($size));
#    print "going to $random_pos\n";
    seek( $file1, $random_pos, 0);
    while ( <$file1> ) {
      last if ( $_ =~ /^\@\w+:\d+/);
    }
    
    if ($_ &&  $_ =~ /^\@\w+:\d+/ && !$used{ $_ }) {
      my $name = $_;
      my $seq  = <$file1>;
      my $str  = <$file1>;
      my $qual = <$file1>;
            
      print $out1 join("", $name, $seq, $str, $qual);
      $read += length("$name$seq$str$qual");
      $used{$infile1}++;
      $name1 = $name;
    }

    if ( $infile2 && $name1) {
      
      $name1 =~ s/\/\d\n//;
      
      my %seen;
      
    REREAD:
      
      seek( $file2, $random_pos, 0);
      while ( <$file2> ) {
	last if ( $_ =~ /^\@\w+:\d+/);
      }

      if ($_ &&  $_ =~ /^\@\w+:\d+/ && !$used{ $_ }) {
	my $name = $_;
	$name2 = $name;
	$name2 =~ s/\/\d\n//;

	next if ($seen{$name2} && $seen{$name2} > 10);


	if ( $name1 eq $name2 ) {
#	  print "$name1 eq $name2\n";
	  my $seq  = <$file2>;
	  my $str  = <$file2>;
	  my $qual = <$file2>;
	  
	
	  print $out2 join("", $name, $seq, $str, $qual);
	  $read += length("$name$seq$str$qual");
	}
	elsif ($name1 lt $name2 ) {
#	  print "$name1 lt $name2\n";
	  $random_pos -= 100;
	  $seen{ $name2 }++;
	  goto REREAD;
	}
	else {
#	  print "$name1 gt $name2\n";
	  $random_pos = tell($file2);
	  $seen{ $name2 }++;
	  goto REREAD;
	}
      }
    }
  }
}



# 
# 
# 
# Kim Brugger (18 Aug 2010)
sub usage {
  
  $0 =~ s/.*\///;
  print "There are three running modes, one where the seq is sampled and aligned and the other from a pre aligned file, finally just qv and base dist.\n";
  print "USAGE RAW: $0 -1 [first fq file] -2 [second fq] -R[eference bwa formatted] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  print "USAGE RAW: $0 -b[am file] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  print "USAGE RAW: $0 -1 [first fq file] -2 [second fq] -n[o mapping] -s[ample size] -o[ut file, otherwise a default will be made] -p[latform]\n";
  exit;

}
