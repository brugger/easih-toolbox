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
getopts('f:b:o:hmu', \%opts);

my $fq_file       = $opts{f};
my $bam_file      = $opts{b};
my $out_file      = $opts{o} || 'test';
my $mapped_only   = $opts{m} || 0;
my $unmapped_only = $opts{u} || 0;
my $combined      = $opts{c} || 0;


my $samtools  = '/usr/local/bin/samtools';

my ( %base_qual_dist, @base_dist, @base_qual, $reads, $Q_count, $mapped, $unmapped);
my ( %mapped_base_qual_dist, @mapped_base_dist, @mapped_base_qual, $mapped_reads, $mapped_Q_count);
my ( %unmapped_base_qual_dist, @unmapped_base_dist, @unmapped_base_qual, $unmapped_reads, $unmapped_Q_count);

my $Q_min = 20;

read_fq($fq_file) if ( $fq_file );
read_bam($bam_file) if ( $bam_file );


report();




# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub report {

  my $infile = $bam_file || $fq_file;

  my ($out, $R);

  if (0) {
  
  open ( $out, "> $out_file.MeanQual ") || die "Could not open '$out_file.MeanQual': $!\n";
  for(my $i = 0; $i < @base_qual; $i++ ) {
    printf $out ("%d\t%.2f\n", $i+1, $base_qual[ $i ] / $reads);
  }
  close( $out );

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R  "meanQual = read.table('$out_file.MeanQual')\n";
  print $R "png('$out_file.MeanQual.png')\n";
  print $R "plot(meanQual, main='$infile: Quality By Cycle', type='h', xlab='Cycle', ylab='Mean Quality')\n";
  print $R "abline(h=15, col='red')\n";
  print $R "dev.off()\n";
  close ($R);
      

  exit;


  open ($out, "> $out_file.QualHist ") || die "Could not open '$out_file.QualHist': $!\n";
  foreach my $value ( sort { $a <=> $b } keys %base_qual_dist ) {
    print  $out "$value\t$base_qual_dist{ $value }\n";
  }
  close( $out );

  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R  "qualHist = read.table('$out_file.QualHist')\n";
  print $R "png('$out_file.QualHist.png')\n";
  print $R "plot(qualHist, main='$infile:  Quality Score Distribution', type='h', xlab='Quality score', ylab='Observations')\n";
  print $R "dev.off()\n";
  close ($R);


  exit;

  }
  

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

  print $out "\t"  . join("\t", @pos) . "\n";
  print $out "A\t" . join("\t", @As) . "\n";
  print $out "C\t" . join("\t", @Cs) . "\n";
  print $out "G\t" . join("\t", @Gs) . "\n";
  print $out "T\t" . join("\t", @Ts) . "\n";
  print $out "N\t" . join("\t", @Ns) . "\n";

  close ($out);



  open ( $R, " | R --vanilla --slave ") || die "Could not open R: $!\n";
  print $R "png('$out_file.BaseDist.png')\n";
  print $R "baseDist = read.table('$out_file.BaseDist')\n";
  print $R "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
  print $R "barplot(as.matrix(baseDist), col=c('green', 'blue', 'black','red', 'grey'), axis.lty=1, main='$infile: Base distribution', xlab='Cycle', ylab='Distribution')\n";
  print $R "legend(45, 90, rev(rownames(baseDist)), fill=rev(c('green', 'blue', 'black','red', 'grey')))\n";
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

    if ( $flags & 0x0004

    analyse( $sequence, $quality, $insert_size);
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
    
    
  }
  close( $fastq);


}



# 
# 
# 
# Kim Brugger (16 Jul 2010)
sub analyse {
  my ( $seq, $qual, $isize) = @_;

  $reads++;
  
  my @seq = split("", $seq);
  my @qual = split("", $qual);

  for(my $i = 0; $i < @seq; $i++) {
    my $base = uc($seq[$i]);
    my $base_qual = ord($qual[$i]) - 33;


    $base_dist[$i]{$base}++;

    $base_qual[$i] += $base_qual;
    $base_qual_dist{$base_qual}++;

    $Q_count++ if ( $Q_min <= $base_qual );
  }

}




