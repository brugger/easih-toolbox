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

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
my $DYNAMIC_LIB_PATHS = 1;
BEGIN {
  if ( $DYNAMIC_LIB_PATHS ) {
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    use lib '/home/kb468/easih-toolbox/modules/';
  }
}

use EASIH::QC;
use EASIH::Misc;

my %opts;
getopts('b:c:f:q:hs:p:r', \%opts);

my $fastq_file    = $opts{f};
my $csfasta_file  = $opts{c};
my $qual_file     = $opts{q};
my $bam_file      = $opts{b};
my $sample_size   = $opts{s} || 20;
my $platform      = uc($opts{p}) || usage();
my $random_sample = $opts{r} || 0;
#sample limit MB
EASIH::QC::sample_size($sample_size);
EASIH::QC::random_sample(1) if ( $random_sample );


my ( $tmp_dir, $tmp_file) = EASIH::Misc::tmp_dir_file();

usage() if ( ($fastq_file   && $csfasta_file) || 
	     ($fastq_file   && $qual_file) || 
	     ($fastq_file   && $bam_file) || 
#	     ($csfasta_file && $qual_file) || 
	     ($csfasta_file && $bam_file) || 
	     ($qual_file    && $bam_file) || 
	     ! $platform);


$platform = "TORRENT" if ( $platform eq "ION");

my %platforms = ( "ILLUMINA" => "Illumina GA2X",
		  "SOLID"    => "SOLiD v4",
		  "TORRENT"  => "Ion Torrent");
# 
# 
# Kim Brugger (03 Feb 2011)
sub usage {
  $0 =~ s/.*\///;
  print "USAGE: $0 -f[astq file] -c[sfasta file] -q[ual file] -b[am file] -p[latform, either SOLID or ILLUMINA or TORRENT] -s<ample size (in Mbases)> -r<andom sampling>\n";
  exit -1;
}



#print "tmp stuff :: $tmp_dir, $tmp_file \n";

my $pdflatex = EASIH::Misc::find_program('pdflatex');
die "pdflatex is not installed on this system\n" if ( !$pdflatex);

my $cwd      = `pwd`;
chomp($cwd);


my ($QC, $base_name, $infile);

if ( $fastq_file || $csfasta_file || $qual_file ) {

  if ( $fastq_file ) {
    ($QC, undef) = EASIH::QC::fastQC( $fastq_file );
    $base_name = $fastq_file;
    $infile    = $fastq_file;
    if ( $base_name =~ /\.\// || $base_name !~ /^\//) {
      $base_name =~ s/\.\///;
      my $cwd = `pwd`;
      chomp $cwd;
      $base_name = "$cwd/$base_name";
    }
    $base_name =~ s/(.*?\.[12])\..*/$1/ ||  $base_name =~ s/(.*?)\..*/$1/; 
    $base_name = "$cwd/$base_name" if ( $base_name !~ /\//);
    EASIH::QC::make_plots($QC, "$tmp_dir/$tmp_file", $base_name);
  }

  if ( $csfasta_file ) {
    $QC = EASIH::QC::csfastaQC( $csfasta_file);
    $base_name = $csfasta_file;
    $infile    = $csfasta_file;
    if ( $base_name =~ /\.\// || $base_name !~ /^\//) {
      $base_name =~ s/\.\///;
      my $cwd = `pwd`;
      chomp $cwd;
      $base_name = "$cwd/$base_name";
    }
    $base_name =~ s/(.*?)\..*/$1/;
    $base_name = "$cwd/$base_name" if ( $base_name !~ /\//);
    EASIH::QC::make_plots($QC, "$tmp_dir/$tmp_file", $base_name);
  }

  if ( $qual_file ) {
    $QC = EASIH::QC::qualQC( $qual_file );


    $base_name = $qual_file;
    $infile    = $qual_file;
    if ( $base_name =~ /\.\// || $base_name !~ /^\//) {
      $base_name =~ s/\.\///;
      my $cwd = `pwd`;
      chomp $cwd;
      $base_name = "$cwd/$base_name";
    }
    $base_name =~ s/(.*?)\..*/$1/;
    $base_name = "$cwd/$base_name" if ( $base_name !~ /\//);

    EASIH::QC::make_plots($QC, "$tmp_dir/$tmp_file", $base_name);
  }


  open (my $out, ">  $tmp_dir/$tmp_file.tex") || die "Could not open 'outfile': $!\n";
  print $out latex_header();
  print $out latex_summary($infile, uc($platform), $sample_size, $$QC{reads} || $$QC{quals}, $$QC{Q30}, $$QC{perc_dup}, $$QC{mappability}, $$QC{ACsplit}, $$QC{partial_adaptor} );
  print $out latex_QV("$tmp_dir/$tmp_file\_BaseQual.pdf", "$tmp_dir/$tmp_file\_QualHist.pdf");
  print $out latex_dups("$tmp_dir/$tmp_file\_DupHist.pdf", $$QC{duplicates}) if ( $platform eq "ILLUMINA" || $platform eq "TORRENT");
  print $out latex_pam("$tmp_dir/$tmp_file\_PAM.pdf", $$QC{partial_adaptor}) if ( $platform eq "ILLUMINA" || $platform eq "TORRENT");
  print $out latex_GC("$tmp_dir/$tmp_file\_BaseDist.pdf", "$tmp_dir/$tmp_file\_GC.pdf");
  print $out latex_tail();

#  $base_name = '/home/kb468/projects/QC/tyt';

  make_pdf($tmp_dir, "$tmp_file.tex", "$base_name.pdf");
  system "rm -rf $tmp_dir";

#  print "rm -rf $tmp_dir $platform\n";

#  print "Report: $base_name.pdf\n";
}
elsif ( $bam_file ) {
  my ($QC1, $QC2) = EASIH::QC::bamQC( $bam_file );

  $base_name = $bam_file;
  $infile    = $bam_file;
  $base_name =~ s/(.*?)\..*/$1/;
  $base_name = "$cwd/$base_name" if ( $base_name !~ /\//);
  
  EASIH::QC::make_plots($QC1, "$tmp_dir/$tmp_file", $base_name);

  open (my $out, ">  $tmp_dir/$tmp_file.tex") || die "Could not open outfile: $!\n";
  print $out latex_header();
  print $out latex_summary($infile, uc($platform), $sample_size, $$QC1{reads} || $$QC1{quals}, $$QC1{Q30}, $$QC1{perc_dup}, $$QC1{mappability}, $$QC1{ACsplit}, $$QC{partial_adaptor} );
  print $out latex_QV("$tmp_dir/$tmp_file\_BaseQual.pdf", "$tmp_dir/$tmp_file\_QualHist.pdf");
  print $out latex_dups("$tmp_dir/$tmp_file\_DupHist.pdf", $$QC1{duplicates});
  print $out latex_pam("$tmp_dir/$tmp_file\_PAM.pdf", $$QC{partial_adaptor});
  print $out latex_GC("$tmp_dir/$tmp_file\_BaseDist.pdf", "$tmp_dir/$tmp_file\_GC.pdf");
  print $out latex_tail();
  make_pdf($tmp_dir, "$tmp_file.tex", "$base_name.1.pdf");
  system "rm -rf $tmp_dir/*";
  print "Report: $base_name.1.pdf\n";

  EASIH::QC::make_plots($QC2, "$tmp_dir/$tmp_file", $base_name);
  open ( $out, ">  $tmp_dir/$tmp_file.tex") || die "Could not open outfile: $!\n";
  print $out latex_header();
  print $out latex_summary($infile, uc($platform), $sample_size, $$QC2{reads} || $$QC2{quals}, $$QC2{Q30}, $$QC2{perc_dup}, $$QC2{mappability}, $$QC2{ACsplit}, $$QC{partial_adaptor} );
  print $out latex_QV("$tmp_dir/$tmp_file\_BaseQual.pdf", "$tmp_dir/$tmp_file\_QualHist.pdf");
  print $out latex_dups("$tmp_dir/$tmp_file\_DupHist.pdf", $$QC2{duplicates});
  print $out latex_pam("$tmp_dir/$tmp_file\_PAM.pdf", $$QC{partial_adaptor});
  print $out latex_GC("$tmp_dir/$tmp_file\_BaseDist.pdf", "$tmp_dir/$tmp_file\_GC.pdf");
  print $out latex_tail();
  make_pdf($tmp_dir, "$tmp_file.tex", "$base_name.2.pdf");
  system "rm -rf $tmp_dir";
  print "Report: $base_name.2.pdf\n";
}




# 
# 
# 
# Kim Brugger (03 Feb 2011)
sub make_pdf {
  my ($dir, $tex, $report) = @_;

  system "cd $tmp_dir/; $pdflatex $tex >/dev/null ";
  $tex =~ s/tex/pdf/;
  system "cd $tmp_dir/; mv $tex $report";

}

# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_summary {
  my ($input, $platform, $size, $reads, $Q30, $dups, $mappable, $AC, $partial_adaptor  ) = @_;


  $input =~ s/.*\///;
  
  my $s = q(\section*{Summary}) . "\n";
  
  $s .= q(\begin{table}[!h])."\n";
  $s .= q(\begin{tabular}{|l|l|}\hline)."\n";
  $s .= q(\rowcolor[gray]{.8} Input file & \verb|).$input. q(|\\\\\hline)."\n" if ( $input );
  $s .= q(Platform & \verb|).$platforms{$platform}. q(|\\\\\hline)."\n" if ( $input );
  $s .= q(\rowcolor[gray]{.8} Sample size & ).$size. q( Mbases (Random sampling)\\\\\hline)."\n" if ( $size && $random_sample);
  $s .= q(\rowcolor[gray]{.8} Sample size & ).$size. q( Mbases\\\\\hline)."\n" if ( $size && !$random_sample);
  $s .= q(Nr of reads & ).$reads.q( reads\\\\\hline)."\n" if ( $reads );

  if ( $platform eq "SOLID" ) {

    $s .= latex_coloured_row("green",  'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 > 25);
    $s .= latex_coloured_row("yellow", 'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 <= 25 && $Q30 >= 15);
    $s .= latex_coloured_row("red",    'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 < 15);
    
    $s .= latex_coloured_row("green", "Mappable prediction", $mappable.'\%')  if ( $mappable && $mappable > 40);
    $s .= latex_coloured_row("yellow", "Mappable prediction", $mappable.'\%') if ( $mappable && $mappable <= 40 && $mappable >= 20);
    $s .= latex_coloured_row("red", "Mappable prediction", $mappable.'\%')    if ( $mappable && $mappable < 20);
  }
  elsif ( $platform eq "ILLUMINA" ) {
    $s .= latex_coloured_row("green",  'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 > 90);
    $s .= latex_coloured_row("yellow", 'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 <= 90 && $Q30 >= 70);
    $s .= latex_coloured_row("red",    'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 < 70);

    $s .= latex_coloured_row("green", "Mappable prediction", $mappable.'\%')  if ( $mappable && $mappable > 95);
    $s .= latex_coloured_row("yellow", "Mappable prediction", $mappable.'\%') if ( $mappable && $mappable >= 70 && $mappable <= 95);
    $s .= latex_coloured_row("red", "Mappable prediction", $mappable.'\%')    if ( $mappable && $mappable < 70);

    $s .= latex_coloured_row("green", "Duplicated sequences", $dups.'\%')  if ( $dups && $dups < 1 );
    $s .= latex_coloured_row("yellow", "Duplicated sequences", $dups.'\%') if ( $dups && $dups >= 1 && $dups <= 10 );
    $s .= latex_coloured_row("red", "Duplicated sequences", $dups.'\%')    if ( $dups && $dups > 10 );

    if ( $partial_adaptor ) {
      $s .= latex_coloured_row("green", 'Partial adaptors ', $partial_adaptor.'\%')  if ( $partial_adaptor < 1 );
      $s .= latex_coloured_row("yellow", 'Partial adaptors ', $partial_adaptor.'\%') if ( $partial_adaptor >= 1 && $partial_adaptor <= 10 );
      $s .= latex_coloured_row("red", 'Partial adaptors ', $partial_adaptor.'\%')    if ( $partial_adaptor > 10 );
    }
  }
  elsif ( $platform eq "TORRENT" ) {
    $s .= latex_coloured_row("green",  'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 > 20);
    $s .= latex_coloured_row("yellow", 'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 <= 20 && $Q30 >= 15);
    $s .= latex_coloured_row("red",    'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 < 10);

    $s .= latex_coloured_row("green", "Mappable prediction", $mappable.'\%')  if ( $mappable && $mappable > 95);
    $s .= latex_coloured_row("yellow", "Mappable prediction", $mappable.'\%') if ( $mappable && $mappable >= 70 && $mappable <= 95);
    $s .= latex_coloured_row("red", "Mappable prediction", $mappable.'\%')    if ( $mappable && $mappable < 70);

    $s .= latex_coloured_row("green", "Duplicated sequences", $dups.'\%')  if ( $dups && $dups < 1 );
    $s .= latex_coloured_row("yellow", "Duplicated sequences", $dups.'\%') if ( $dups && $dups >= 1 && $dups <= 10 );
    $s .= latex_coloured_row("red", "Duplicated sequences", $dups.'\%')    if ( $dups && $dups > 10 );

    if ( $partial_adaptor ) {
      $s .= latex_coloured_row("green", 'Partial adaptors ', $partial_adaptor.'\%')  if ( $partial_adaptor < 1 );
      $s .= latex_coloured_row("yellow", 'Partial adaptors ', $partial_adaptor.'\%') if ( $partial_adaptor >= 1 && $partial_adaptor <= 10 );
      $s .= latex_coloured_row("red", 'Partial adaptors ', $partial_adaptor.'\%')    if ( $partial_adaptor > 10 );
    }
  }


  if ( $AC ) {
    $s .= latex_coloured_row("green", 'Avg \% AC ', $AC.'\%')  if ( $AC >= 45 && $AC <= 55 );
    $s .= latex_coloured_row("yellow", 'Avg \% AC ', $AC.'\%') if ( $AC < 45 && $AC > 40 );
    $s .= latex_coloured_row("yellow", 'Avg \% AC ', $AC.'\%') if ( $AC < 60 && $AC > 55 );
    $s .= latex_coloured_row("red", 'Avg \% AC ', $AC.'\%')    if ( $AC < 40 || $AC > 60 );
  }



  $s .= q(\end{tabular}) . "\n";
  $s .= q(\end{table}) . "\n";

  $s .= "\n\n\n";
  $s .= q(\begin{table}[!h])."\n";
  $s .= q(\begin{tabular}{|l|>{\columncolor{red}}c|>{\columncolor{yellow}}c|>{\columncolor{green}}c|}\hline)."\n";
  $s .= q(Bases $>=$Q30         & $<15 \%$ &  $15-25 \%$ & $>25 \%$\\\\\hline )."\n" if ( $Q30 && $platform eq "SOLID");
  $s .= q(Bases $>=$Q30      & $<70 \%$ &  $70-90 \%$ & $>90 \%$\\\\\hline )."\n" if ( $Q30 && $platform eq "ILLUMINA");
  $s .= q(Mappable prediction     & $<20 \%$ &  $20-40 \%$ & $>40 \%$\\\\\hline)."\n" if ($mappable && $platform eq "SOLID");
  $s .= q(Mappable prediction  & $<70 \%$ &  $70-95 \%$ & $>95 \%$\\\\\hline)."\n" if ($mappable && $platform eq "ILLUMINA");
  $s .= q(Duplicates sequences & $>10 \%$ &  $1-10  \%$ & $<1 \%$\\\\\hline)."\n" if ( $dups);
  $s .= q(Avg \% AC & $<40\%$ or $>60\%$ & $40-45\%$ or $55-60\%$ & $45-55 \%$\\\\\hline)."\n" if ( $AC );
  $s .= q(Partial adaptors & $>10 \%$ &  $1-10  \%$ & $<1 \%$\\\\\hline)."\n" if ( $partial_adaptor );

  $s .= q(\end{tabular})."\n";
  $s .= q(\end{table})."\n";

  $s .= q(\newpage)."\n";

  return $s;

}



# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_coloured_row {
  my ($colour, @values) = @_;

  return '\rowcolor{' . $colour. '}' . join(" & ", @values) . ' \\\\\hline'."\n";
}



# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_pam {
  my ($PAM, $partial_adaptor) = @_;

  print "$PAM\n";

  return if ( !$PAM || ! -e $PAM );

  my $s = q|\section*{Partial mapped adaptors}|."\n";

  $s  .= q(\begin{figure}[!h]) . " \n";
  $s  .= q(\centering) . " \n";
  $s  .= q(\includegraphics[scale=0.5]{). $PAM .q(}) . " \n";
  $s  .= q(\end{figure}) . " \n";

  $s  .= q(The figure shows the number of hits to adaptor sequences pr base (bars on the X-axis) with the number of) . " \n";
  $s  .= q(observations (Y-axis). ) . " \n";

  $s  .= q(\newpage) . " \n";

  return $s;

}



# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_dups {
  my ($DupHist, $duplicates) = @_;

  return if ( !$DupHist && ! $duplicates);
  return if ( ! -e $DupHist );

  my $s = q|\section*{Duplications}|."\n";

  $s  .= q(\begin{figure}[!h]) . " \n";
  $s  .= q(\centering) . " \n";
  $s  .= q(\includegraphics[scale=0.5]{). $DupHist .q(}) . " \n";
  $s  .= q(\end{figure}) . " \n";

  $s  .= q(The figure shows the number of sequences (bars on the X-axis) with the number of) . " \n";
  $s  .= q(observations (Y-axis). A further breakdown of the sequences, and their) . " \n";
  $s  .= q(potential origin are shown on the next page.) . " \n";

  $s  .= q(\newpage) . " \n";

  $s  .= q(\landscape) . " \n";

  $s  .= q(\begin{table}[!h]) . " \n";
  $s  .= q(\begin{tabular}{|l|l|l|} \\hline) . " \n";
  $s  .= q(\rowcolor[gray]{.8}Sequence & Percentage & Possible source \\\\\hline) . " \n";

  
  my $max = 25;
  foreach my $read ( sort { $$duplicates{$b}{count} <=>  $$duplicates{$a}{count} } keys %$duplicates ) {
    my $type = EASIH::QC::check_for_adaptors($read);
    $s  .= join(" & ", q(\tiny{) . $read. q(}), $$duplicates{$read}{percent}.'\%',q(\tiny{) .$type.q(})). q(  \\\\\hline) . " \n";
    last if ( ! $max--);
  }
  $s  .= q(\end{tabular}) . " \n";
  $s  .= q(\end{table}) . " \n";

  $s  .= q(\endlandscape) . " \n";
  $s  .= q(\newpage)."\n";

  return $s;

}


# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_QV {
  my ($BaseQual, $QualHist) = @_;


  return "" if ( ! $BaseQual && ! $QualHist );
  return if ( ! -e $BaseQual && ! -e $QualHist );

  my $s = q|\section*{Quality Values}| . "\n";


  $s  .= q(The quality values normally range between 0 and 40. These are phred
like scores. The Phred table shows how the scores relates to the expected
number of error, for more information see:
\verb|http://en.wikipedia.org/wiki/Phred_quality_score|. As a rule of
thumb the mean QV should be about 25-30.) . " \n";


  $s .= q(\begin{table}[!h])."\n";
  $s .= q(\begin{tabular}{|l|l|l|}\hline)."\n";

  $s .= q(\rowcolor[gray]{.8}Phred Quality Score & Probability of incorrect base call &	Base call accuracy \\\\\hline)."\n";
  $s .= q(10 & 1 in 10 & 90 \% \\\\\hline)."\n";
  $s .= q(20 & 1 in 100 & 99 \% \\\\\hline)."\n";
  $s .= q(30 & 1 in 1000 & 99.9 \% \\\\\hline)."\n";
  $s .= q(40 & 1 in 10000 & 99.99 \% \\\\\hline)."\n";
  $s .= q(50 & 1 in 100000 & 99.999 \% \\\\\hline)."\n";
  $s .= q(\end{tabular})."\n";
  $s .= q(\end{table})."\n";

  $s .= q(\begin{figure}[ht])."\n";
  $s .= q(\hspace{-0.5cm})."\n";
  $s .= q(\begin{minipage}[b]{0.54\linewidth})."\n";
  $s .= q(\centering)."\n";
  $s .= q(\includegraphics[scale=0.5]{). $BaseQual .q(})."\n" if ($BaseQual && -e $BaseQual);
  $s .= q(\end{minipage})."\n";
  $s .= q(\hspace{0.5cm})."\n";
  $s .= q(\begin{minipage}[b]{0.54\linewidth})."\n";
  $s .= q(\centering)."\n";
  $s .= q(\includegraphics[scale=0.5]{). $QualHist .q(})."\n" if ($QualHist && -e $QualHist);
  $s .= q(\end{minipage})."\n";
  $s .= q(\end{figure})."\n";

  $s .= q(\newpage)."\n";
  
  return $s;
}




# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_GC {
  my ($BaseDist, $GC) = @_;

  return if ( !$BaseDist && ! $GC);
  return if ( ! -e $BaseDist && ! -e $GC);
  

  my $s = q|\section*{Base/colour distribution and GC\%}| ."\n";


  $s .= q(\begin{figure}[ht])."\n";
  $s .= q(\hspace{-0.5cm})."\n";
  $s .= q(\begin{minipage}[b]{0.54\linewidth})."\n";
  $s .= q(\centering)."\n";
  $s .= q(\includegraphics[scale=0.5]{).$BaseDist.q(})."\n" if ( $BaseDist && -e $BaseDist );
  $s .= q(\end{minipage})."\n";
  $s .= q(\hspace{0.5cm})."\n";
  $s .= q(\begin{minipage}[b]{0.54\linewidth})."\n";
  $s .= q(\centering)."\n";
  $s .= q(\includegraphics[scale=0.5]{).$GC.q(})."\n" if ( $GC && -e $GC );
  $s .= q(\end{minipage})."\n";
  $s .= q(\end{figure})."\n";

  $s .= q(The base/colour distribution plot shows the percentages of each
base/colour at a given sequencing cycle. For samples with high
variation the splits should \~25\% for each colour show, if the
samples are low variation samples like amplicons, CHIP-seq this is
normally not the case. The base/colour shown corresponds to A/0:
green, C/1: yellow, G/2: black, T/3: blue, N/.: grey.

The GC distribution should from a high variation sample display a bell
shaped curve. Low variation samples like Amplicons, CHIP-seq will most
likely not be this. Odd bumps and shoulders on the curve could
originate from contaminating or sequencing adaptors.

\end{document}

);

}



# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_header {


  return q(
\include{mixture}

\setlength{\emergencystretch}{2em}

\documentclass[12pt, a4paper, oneside]{report}
%%\special{papersize=210mm,2970mm}
\usepackage[hmargin=3cm,vmargin=3.5cm]{geometry}

\linespread{1.0} %1.3 == 1.5,  1.6 == double
\usepackage{graphicx}
\usepackage{fancyhdr}
\pagestyle{fancy}

\usepackage{colortbl}
\usepackage{color}
\usepackage{array}
\usepackage{lscape}

\fancyhf{}
\fancypagestyle{plain}{%
\fancyhead{} % get rid of headers on plain pages
\renewcommand{\headrulewidth}{0pt} % and the line
}
%%\lhead{} QC report}
\rhead{EASIH QC pipeline}
\begin{document}

)

}



# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_tail {
  return q(
\end{document}
);
}



__END__

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
