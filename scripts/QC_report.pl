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
my $sample_size   = $opts{s} || 1;

#sample limit MB
EASIH::QC::sample_size($sample_size);

if ( $fastq_file ) {
  my $QC = EASIH::QC::fastQC( $fastq_file, 1 );
  $fastq_file =~ s/(.*?.[12])\..*/$1/ ||  $fastq_file =~ s/(.*?)\..*/$1/; 
  EASIH::QC::make_plots($QC, $fastq_file);
  print "Predicted mappability : $$QC{mappability} %\n" if ( $$QC{mappability} );

  open (my $out, "> /home/kb468/projects/QC/latex/test.tex") || die "Could not open outfile: $!\n";
 print $out latex_header($fastq_file);
  print $out latex_summary($opts{f}, "ILLUMINA", $sample_size, $$QC{reads}, $$QC{Q30}, $$QC{perc_dup}, $$QC{mappability}, $$QC{ACsplit} );
  print $out latex_QV("/data/A06/raw/A0600028_1_BaseQual.pdf", "/data/A06/raw/A0600028_1_BaseQualDist_QualHist.pdf");
  print $out latex_dups("/data/A06/raw/A0600028_1_Duplicates_DupHist.pdf", $$QC{duplicates});
  print $out latex_GC("/data/A06/raw/A0600028_1_BaseDist.pdf", "/data/A06/raw/A0600028_1_GC.pdf");
  print $out latex_tail();

}

if ( $csfastq_file ) {
  my $QC = EASIH::QC::csfastaQC( $csfastq_file);
  $csfastq_file =~ s/(.*?)\..*/$1/;
  EASIH::QC::make_plots($QC, $csfastq_file);
}

if ( $qual_file ) {
  my $QC = EASIH::QC::qualQC( $qual_file );
  $qual_file =~ s/(.*?)\..*/$1/;
  EASIH::QC::make_plots($QC, $qual_file);
  print "Predicted mappability : $$QC{mappability} %\n" if ( $$QC{mappability} );
  print "+Q30: $$QC{Q30}\n" if ($$QC{Q30});

  open (my $out, "> /home/kb468/projects/QC/latex/test.tex") || die "Could not open outfile: $!\n";

  print $out latex_header($qual_file);
  print $out latex_summary($opts{q}, "SOLID", $sample_size, $$QC{reads}, $$QC{Q30}, undef, $$QC{mappability}, undef );
  print $out latex_QV("$qual_file\_BaseQual.pdf", "$qual_file\_BaseQualDist_QualHist.pdf");
#  print $out latex_dups($DupHist, $dup_reads);


  print $out latex_tail();

}

if ( $bam_file ) {
  my ($QC1, $QC2) = EASIH::QC::bamQC( $bam_file );
  $bam_file =~ s/(.*?)\..*/$1/;
  
  EASIH::QC::make_plots($QC1, "$bam_file.1");
  print "Predicted mappability : $$QC1{mappability} %\n" if ( $$QC1{mappability} );

  EASIH::QC::make_plots($QC2, "$bam_file.2");
  print "Predicted mappability : $$QC2{mappability} %\n" if ( $$QC2{mappability} );
}



# 
# 
# 
# Kim Brugger (02 Feb 2011)
sub latex_summary {
  my ($input, $platform, $size, $reads, $Q30, $dups, $mappable, $AC,  ) = @_;

  $input =~ s/.*\///;
  
  my $s = q(\section*{Summary}) . "\n";


  $s .= q(\begin{table}[!h])."\n";
  $s .= q(\begin{tabular}{|l|l|}\hline)."\n";
  $s .= q(\rowcolor[gray]{.8} Input file & \verb|).$input. q(|\\\\\hline)."\n" if ( $input );
  $s .= q(Sample size & ).$size. q( Mbases\\\\\hline)."\n" if ( $size);
  $s .= q(\rowcolor[gray]{.8}Nr of reads & ).$reads.q( reads\\\\\hline)."\n" if ( $reads );

  if ( $platform eq "SOLID" ) {

    $s .= latex_coloured_row("green",  'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 > 40);
    $s .= latex_coloured_row("yellow", 'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 <= 40 && $Q30 >= 20);
    $s .= latex_coloured_row("red",    'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 < 20);
    
    $s .= latex_coloured_row("green", "Mappable prediction", $mappable.'\%') if ( $mappable && $mappable > 40);
    $s .= latex_coloured_row("yellow", "Mappable prediction", $mappable.'\%') if ( $mappable && $mappable <= 40 && $mappable >= 20);
    $s .= latex_coloured_row("red", "Mappable prediction", $mappable.'\%')   if ( $mappable && $mappable < 20);
  }
  elsif ( $platform eq "ILLUMINA" ) {
    $s .= latex_coloured_row("green",  'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 > 90);
    $s .= latex_coloured_row("yellow", 'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 <= 90 && $Q30 >= 70);
    $s .= latex_coloured_row("red",    'Bases $>=$Q30', $Q30.'\%') if ( $Q30 && $Q30 < 70);

    $s .= latex_coloured_row("green", "Mappable prediction", $mappable.'\%') if ( $mappable && $mappable > 95);
    $s .= latex_coloured_row("yellow", "Mappable prediction", $mappable.'\%') if ( $mappable && $mappable <= 70 && $mappable >= 95);
    $s .= latex_coloured_row("red", "Mappable prediction", $mappable.'\%')   if ( $mappable && $mappable < 70);
  }

  $s .= latex_coloured_row("green", "Duplicated sequences", $dups.'\%') if ( $dups && $dups < 1 );
  $s .= latex_coloured_row("yellow", "Duplicated sequences", $dups.'\%') if ( $dups && $dups >= 1 && $dups <= 1 );
  $s .= latex_coloured_row("red", "Duplicated sequences", $dups.'\%') if ( $dups && $dups > 10 );

  $s .= latex_coloured_row("green", 'Avg \% AC ', $AC.'\%') if ( $AC && $AC >= 45 && $AC <= 55);
  $s .= latex_coloured_row("yellow", 'Avg \% AC ', $AC.'\%') if ( $AC && $AC < 45 && $AC > 40);
  $s .= latex_coloured_row("yellow", 'Avg \% AC ', $AC.'\%') if ( $AC && $AC < 60 && $AC > 55);
  $s .= latex_coloured_row("red", 'Avg \% AC ', $AC.'\%') if ( $AC && $AC < 40 || $AC > 60 );

  $s .= q(\end{tabular}) . "\n";
  $s .= q(\end{table}) . "\n";

  $s .= "\n\n\n";
  $s .= q(\begin{table}[!h])."\n";
  $s .= q(\begin{tabular}{|l|>{\columncolor{red}}c|>{\columncolor{yellow}}c|>{\columncolor{green}}c|}\hline)."\n";
  $s .= q(Bases $>=$Q30 (SOLiD)        & $<20 \%$ &  $20-40 \%$ & $>40 \%$\\\\\hline )."\n" if ( $Q30);
  $s .= q(Bases $>=$Q30 (Illumina)     & $<70 \%$ &  $70-90 \%$ & $>90 \%$\\\\\hline )."\n" if ( $Q30);
  $s .= q(Mappable prediction (SOLiD)    & $<20 \%$ &  $20-40 \%$ & $>40 \%$\\\\\hline)."\n" if ($mappable && $platform eq "SOLID");
  $s .= q(Mappable prediction (Illumina) & $<70 \%$ &  $70-95 \%$ & $>95 \%$\\\\\hline)."\n" if ($mappable && $platform eq "ILLUMINA");
  $s .= q(Duplicates sequences & $>10 \%$ &  $1-10  \%$ & $<1 \%$\\\\\hline)."\n" if ( $dups);
  $s .= q(Avg \% AC & $<40\%$ or $>60\%$ & $40-45\%$ or $55-60\%$ & $45-55 \%$\\\\\hline)."\n" if ( $AC );

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
sub latex_dups {
  my ($DupHist, $dup_reads) = @_;

  return if ( !$DupHist && ! $dup_reads);

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
  $s  .= q(\begin{tabular}{|l|l|l|} \hline) . " \n";
  $s  .= q(\rowcolor[gray]{.8}Sequence & Percentage & Possible source \\\hline) . " \n";

#  foreach my $read ( @$dup_reads ) {
#    $s  .= join(" & ", q(\tiny{) . $$read{seq}. q(}), $$read{perc}.'\%',$$read{type}). q(  \\\hline) . " \n";
#  }
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
  $s .= q(\includegraphics[scale=0.5]{). $BaseQual .q(})."\n" if ($BaseQual);
  $s .= q(\end{minipage})."\n";
  $s .= q(\hspace{0.5cm})."\n";
  $s .= q(\begin{minipage}[b]{0.54\linewidth})."\n";
  $s .= q(\centering)."\n";
  $s .= q(\includegraphics[scale=0.5]{). $QualHist .q(})."\n" if ($QualHist);
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

  my $s = q|\section*{Base/colour distribution and GC\%}| ."\n";


  $s .= q(\begin{figure}[ht])."\n";
  $s .= q(\hspace{-0.5cm})."\n";
  $s .= q(\begin{minipage}[b]{0.54\linewidth})."\n";
  $s .= q(\centering)."\n";
  $s .= q(\includegraphics[scale=0.5]{).$BaseDist.q(})."\n";
  $s .= q(\end{minipage})."\n";
  $s .= q(\hspace{0.5cm})."\n";
  $s .= q(\begin{minipage}[b]{0.54\linewidth})."\n";
  $s .= q(\centering)."\n";
  $s .= q(\includegraphics[scale=0.5]{).$GC.q(})."\n";
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
  my ($name) = @_;

  $name =~ s/.*\///;

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
%%\lhead{).$name.q(\end{verbatim} QC report}
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
