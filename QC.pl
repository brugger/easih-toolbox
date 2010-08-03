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
getopts('f:b:o:hmusp:', \%opts);

my $fq_file       = $opts{f};
my $bam_file      = $opts{b};
my $out_file      = $opts{o};
my $mapped_only   = $opts{m} || 0;
my $prefix        = $opts{p} || "";
my $unmapped_only = $opts{u} || 0;
my $solid         = $opts{s} || 0;

my $limit         = 10000;

my $samtools  = '/usr/local/bin/samtools';

my $MAPPED_READ    =  0;
my $UNMAPPED_READ  =  4;

my ( %base_qual_dist, @base_dist, @base_qual);
my ( %mapped_base_qual_dist, @mapped_base_dist, @mapped_base_qual);
my ( %unmapped_base_qual_dist, @unmapped_base_dist, @unmapped_base_qual);
my ( %GC );

my ($reads, $mapped_reads, $unmapped_reads) = (0,0,0);

my $Q_min = 20;

read_fq($fq_file) if ( $fq_file );
read_bam($bam_file) if ( $bam_file );

die "no input, read the code\n" if ( ! $fq_file && ! $bam_file );



report();

# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub report {

  my $infile = $bam_file || $fq_file;
  $infile =~ s/.*\///;

  $infile = $prefix . $infile;

  $out_file ||= $infile;

  my ($Q_min_count) = (0,0);

  my ($out, $R);
  
  my $max_qual = 0;
  
  open ( $out, "> $out_file.MeanQual ") || die "Could not open '$out_file.MeanQual': $!\n";
  for(my $i = 0; $i < @base_qual; $i++ ) {
    my @values = ($i+1, $base_qual[ $i ] / $reads);
    $max_qual = $base_qual[ $i ] / $reads if ( $max_qual < $base_qual[ $i ] / $reads);

    if ( ! $mapped_reads ) {
      push @values, 0;
    }
    else {
      push @values, ($mapped_base_qual[ $i ] || 0)/ $mapped_reads;
      $max_qual = $mapped_base_qual[ $i ] / $mapped_reads if ( $max_qual < $mapped_base_qual[ $i ] / $mapped_reads);
    }

    if ( ! $unmapped_reads ) {
      push @values, 0;
    }
    else {
      push @values, ($unmapped_base_qual[ $i ] || 0)/ $unmapped_reads;
      $max_qual = $unmapped_base_qual[ $i ] / $reads if ( $max_qual < $unmapped_base_qual[ $i ] / $unmapped_reads);
    }


    print $out join("\t", @values) . "\n";
  }
  close( $out );
  
  $max_qual += 5;

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "meanQual = read.table('$out_file.MeanQual', check.names=FALSE)\n";
#  print $R "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
  print $R "png('$out_file.MeanQual.png')\n";
  print $R "plot(meanQual[,2], main='$infile: Quality By Cycle', type='h', xlab='Cycle', ylab='Mean Quality', ylim=c(0,$max_qual))\n";
  print $R "lines(meanQual[,3], type='l', col='blue')\n";
  print $R "lines(meanQual[,4], type='l', col='red')\n";
  print $R "abline(h=15, col='red')\n";
  print $R "dev.off()\n";
  close ($R);

      
  open ($out, "> $out_file.QualHist ") || die "Could not open '$out_file.QualHist': $!\n";
  my (@values, @combined, @mapped, @unmapped, %saw);
  foreach my $value ( sort { $a <=> $b } grep(!$saw{$_}++,  keys %base_qual_dist, keys %unmapped_base_qual_dist, keys %mapped_base_qual_dist )) {
    push @values, $value;
    push @combined, ($base_qual_dist{ $value } || 0);
    push @mapped,   ($mapped_base_qual_dist{ $value } || 0);
    push @unmapped, ($unmapped_base_qual_dist{ $value } || 0);

    $Q_min_count += $base_qual_dist{ $value } if ( $Q_min <= $value);
  }

  print $out "\t" . join("\t", @values) . "\n";
  print $out "\t" . join("\t", @combined) . "\n";
  print $out "\t" . join("\t", @mapped) . "\n";
  print $out "\t" . join("\t", @unmapped) . "\n";

  close( $out );

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R  "qualHist = read.table('$out_file.QualHist', header=TRUE, check.names=FALSE)\n";
  print $R "png('$out_file.QualHist.png')\n";
#  print $R "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
  print $R "barplot(as.matrix(qualHist), col=c('grey', 'blue', 'red'), axis.lty=1, main='$infile: Quality distribution', xlab='Quality', ylab='Observations', beside=TRUE)\n";
#  print $R "legend(45, 90, c('All', 'Mapped', 'Unmapped', fill=c('grey', 'blue', 'red')))\n";

#  print $R "plot(qualHist, main='$infile:  Quality Score Distribution', type='h', xlab='Quality score', ylab='Observations')\n";
  print $R "dev.off()\n";
  close ($R);
  
  open ( $out, "> $out_file.BaseDist ") || die "Could not open '$out_file.BaseDist': $!\n";

  my (@As, @Cs, @Gs, @Ts, @Ns, @pos);

  for(my $i = 0; $i < @base_dist; $i++ ) {
    my $Ns = $reads - ($base_dist[ $i ]{A} || 0) - ($base_dist[ $i ]{C} || 0) - ($base_dist[ $i ]{G} || 0) - ($base_dist[ $i ]{T} || 0);

    push @pos, $i + 1;
    push @As, ($base_dist[ $i ]{A} || 0)/$reads*100;
    push @Cs, ($base_dist[ $i ]{C} || 0)/$reads*100;
    push @Gs, ($base_dist[ $i ]{G} || 0)/$reads*100;
    push @Ts, ($base_dist[ $i ]{T} || 0)/$reads*100;
    push @Ns, ($Ns || 0)/$reads*100;

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

  open ($out, " > $infile.report ") || die "Could not open '$infile.report': $!\n";

  printf $out ("Density: %.2f\n", $reads/120);
  close($out);

  open ( $out, "> $out_file.GC ") || die "Could not open '$out_file.GC': $!\n";
  foreach my $key ( sort {$a <=> $b} keys %GC ) {
    print $out "$key\t$GC{$key}\n";
  }
  close( $out );

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "png('$out_file.GC.png')\n";
  print $R "GC = read.table('$out_file.GC')\n";
  print $R "a = hist(GC[,1], breaks=10, plot=FALSE)\n";
  print $R "plot(spline(a\$mids, a\$intensities*1000), type='b', xlab='%GC', ylab='fraction')\n";
  print $R "dev.off()\n";
  close ($R);

}




# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub read_bam {
  my ($infile) = @_;

  open (my $st_pipe, "$samtools view $infile  | ") || die "Could not open samtools pipe: $!";

  while(<$st_pipe>) {
    chomp;
    my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");    


    next if ( $flags & 0x0004 && $mapped_only );

    analyse( $sequence, $quality, $insert_size, ($flags & 0x0004));
  }
  close( $st_pipe );
}







# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub read_fq {
  my ( $infile ) = @_;
  
  open (my $fastq, " $infile") || die "Could not open '$infile': $!\n";

  while (<$fastq>) {

    my $name   = $_;
    chomp( $name );
    my $seq    = <$fastq>;
    chomp( $seq );
    my $strand = <$fastq>;
    chomp( $strand );
    my $qual   = <$fastq>;
    chomp( $qual );
    
    analyse($seq, $qual);
    

#    last if (--$limit <=0 );
    
  }
  close( $fastq);


}



# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub analyse {
  my ( $seq, $qual, $isize, $mapped) = @_;

  $reads++;

  $mapped_reads++   if ( defined $mapped && $mapped ==  $MAPPED_READ );
  $unmapped_reads++ if ( defined $mapped &&  $mapped ==  $UNMAPPED_READ );
      

  
  my @seq = split("", $seq);
  my @qual = split("", $qual);

  my ( $GC, $AT) = (0,0);

  for(my $i = 0; $i < @seq; $i++) {
    my $base = uc($seq[$i]);
    my $base_qual = ord($qual[$i]) - 33;

    $base_dist[$i]{$base}++ if (!$solid || defined $mapped && $mapped ==  $MAPPED_READ );

    $GC++ if ( $base eq 'G' || $base eq 'C');
    $AT++ if ( $base eq 'A' || $base eq 'T');

    $base_qual[$i] += $base_qual;
    $base_qual_dist{$base_qual}++;

    if ( defined $mapped && $mapped == $MAPPED_READ ) {
      $mapped_base_dist[$i]{$base}++;
      $mapped_base_qual[$i] += $base_qual;
      $mapped_base_qual_dist{$base_qual}++;
    }
    elsif ( defined $mapped &&  $mapped ==  $UNMAPPED_READ ) {
      $unmapped_base_dist[$i]{$base}++ if ( ! $solid );
      $unmapped_base_qual[$i] += $base_qual;
      $unmapped_base_qual_dist{$base_qual}++;
    }
  }

  $GC{ int( $GC*100/($GC+$AT))}++;

}




