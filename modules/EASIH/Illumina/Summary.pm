package EASIH::Illumina::Summary;
# 
# For pulling data out of Summary.htm files
# 
# 
# Kim Brugger (17 Nov 2011), contact: kim.brugger@easih.ac.uk


use strict;
use warnings;
use Data::Dumper;
use XML::Simple;



# 
# 
# 
# Kim Brugger (17 Nov 2011)
sub readin_summaries {
  my ($indir) = @_;

  opendir(DIR, "$indir");
  my @files = grep(!/^\./, sort readdir(DIR));
  my %res;

  foreach my $infile ( @files ) {
    $infile = "$indir/$infile";
#$infile = '/seqs/hiseq04/111115_SN1078_0071_BC081DACXX/Data/reports/Summary/read1.xml';
    my $config = XMLin($infile);
  
#foreach my $lane ( keys %{$$config{'Lane'}} ) {
    foreach my $lane ( 6 ) {
      $res{$$config{'Read'}}{$lane}{'Phas'}  = $$config{'Lane'}{$lane}{Phasing};
      $res{$$config{'Read'}}{$lane}{'Preph'} = $$config{'Lane'}{$lane}{Prephasing};
      
      $res{$$config{'Read'}}{$lane}{'FirstCycleInt'}   = $$config{'Lane'}{$lane}{FirstCycleIntPF};
      $res{$$config{'Read'}}{$lane}{'FirstCycleIntSD'} = $$config{'Lane'}{$lane}{FirstCycleIntPFSD};
      
      $res{$$config{'Read'}}{$lane}{'20CyclesIntSD'} = $$config{'Lane'}{$lane}{PrcIntensityAfter20CyclesPFSD};
      $res{$$config{'Read'}}{$lane}{'20CyclesIntSD'} = $$config{'Lane'}{$lane}{PrcIntensityAfter20CyclesPFSD};
      
      $res{$$config{'Read'}}{$lane}{'PrcAlign'}   = $$config{'Lane'}{$lane}{PrcAlign};
      $res{$$config{'Read'}}{$lane}{'PrcAlignSD'} = $$config{'Lane'}{$lane}{PrcAlignSD};
      
      $res{$$config{'Read'}}{$lane}{PrcPFClusters}   = $$config{'Lane'}{$lane}{PrcPFClusters};
      $res{$$config{'Read'}}{$lane}{PrcPFClustersSD} = $$config{'Lane'}{$lane}{PrcPFClustersSD};
      
      $res{$$config{'Read'}}{$lane}{'Clusters'}   = sprintf("%d",$$config{'Lane'}{$lane}{ClustersRaw}*$$config{densityRatio});
      $res{$$config{'Read'}}{$lane}{'ClustersSD'} = sprintf("%d",$$config{'Lane'}{$lane}{ClustersRawSD}*$$config{densityRatio});
      
      $res{$$config{'Read'}}{$lane}{'ClustersPF'}   = sprintf("%d",$$config{'Lane'}{$lane}{ClustersPF}*$$config{densityRatio});
      $res{$$config{'Read'}}{$lane}{'ClustersPFSD'} = sprintf("%d",$$config{'Lane'}{$lane}{ClustersPFSD}*$$config{densityRatio});
    }
    
  }

  return \%res;
}


1;



