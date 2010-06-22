#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (28 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/home/kb468/projects/ensembl-variation/modules/';
use lib '/home/kb468/projects/e57/ensembl/modules/';
use lib '/home/kb468/projects/e57/bioperl-live/';

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

my $species     = "human";
my $buffer_size = 500;
my $host        = 'ensembldb.ensembl.org';
my $user        = 'anonymous';

# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);
# get variation adaptors
my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
my $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');
my $sa = $reg->get_adaptor($species, 'core', 'slice');
my $ga = $reg->get_adaptor($species, 'core', 'gene');

my $bed         = shift;
my $from_36     = 0;
my $links       = 0;
my $full_report = 0;
my $html_out    = 0;

my $dbsnp_link    = 'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=';
my $ens_gene_link = 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=';

my ($mapper, $asma, $csa, $cs_from, $cs_to);
if ( $from_36 ) {

  $asma = $reg->get_adaptor($species, 'core', 'assemblymapper');
  $csa  = $reg->get_adaptor($species, 'core', 'coordsystem');

  $cs_from        = 'NCBI36';
  $cs_to          = 'GRCh37';

  $cs_from = $csa->fetch_by_name('chromosome', $cs_from );
  die "Unknown coord system: $cs_from\n" if ( !$cs_from );
  $cs_to   = $csa->fetch_by_name('chromosome', $cs_to);
  die "Unknown coord system: $cs_to\n" if ( !$cs_to );
  $mapper  = $asma->fetch_by_CoordSystems( $cs_from, $cs_to );

}


my $indels = readin_bed( $bed )      if ( $bed);
foreach my $chr ( sort {$a cmp $b}  keys %$indels ) {
  

  foreach my $start ( sort { $a <=> $b} keys %{$$indels{$chr}} ) {

    my $indel = $$indels{$chr}{ $start };
    my $end   = $$indel{end};
    my $position = "$chr:$start-$end";

    print Dumper( $indel );
    
    my $res = indel_effect($chr, $start, $end, "$$indel{variation}/$$indel{type}");

#    print_oneliner( $position, @$res);

  }

  last;
}





# 
# 
# 
# Kim Brugger (28 May 2010)
sub indel_effect {
  my ( $chr, $start, $end, $allele_string ) = @_;

  $chr =~ s/chr//;

  my  $strand = 1;

  if ( $from_36 ) {
    
    my @res = $mapper->map($chr, $start, $end, $strand, $cs_from);
    foreach my $res ( @res ) {
      if ( $res->isa( 'Bio::EnsEMBL::Mapper::Coordinate' )) {
	my $chr_slice = $sa->fetch_by_seq_region_id($res->id);
#	print "$chr, $start, $end --> " . join("\t", $chr_slice->seq_region_name, $res->start, $res->end ) . "\n";
	$start  = $res->start;
	$end    = $res->end;
	$strand = $res->strand;
	last;
      }
    }
  }

  my @res;
  $allele_string = "C/GGGG";
  print "AS :: $allele_string\n";

  my $slice;
  # first try to get a chromosome
  eval { $slice = $sa->fetch_by_region('chromosome', $chr); };
  
  # if failed, try to get any seq region
  if(!defined($slice)) {
    $slice = $sa->fetch_by_region(undef, $chr);
  }
  
  # if failed, die
  if(!defined($slice)) {
    die("ERROR: Could not fetch slice named $chr\n");
  }	

  my @vfs;
  
  # create a new VariationFeature object
  my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
    -start => $start,
    -end => $end,
    -slice => $slice,           # the variation must be attached to a slice
    -allele_string => $allele_string,
    -strand => $strand,
    -map_weight => 1,
    -adaptor => $vfa,           # we must attach a variation feature adaptor
    -variation_name => $chr.'_'.$start.'_'.$allele_string,
  );
  push @vfs, $new_vf; 
  
  $tva->fetch_all_by_VariationFeatures( \@vfs );
  foreach my $vf (@vfs) {

    # find any co-located existing VFs
    my $existing_vf;
    
    if(defined($new_vf->adaptor->db)) {
      my $fs = $new_vf->feature_Slice;
      if($fs->start > $fs->end) {
	($fs->{'start'}, $fs->{'end'}) = ($fs->{'end'}, $fs->{'start'});
      }
      foreach my $existing_vf_obj(@{$new_vf->adaptor->fetch_all_by_Slice($fs)}) {
	$existing_vf = $existing_vf_obj->variation_name
	    if ($existing_vf_obj->seq_region_start == $new_vf->seq_region_start &&
		$existing_vf_obj->seq_region_end   == $new_vf->seq_region_end );
      }
    }
		
    # the get_all_TranscriptVariations here now just retrieves the
    # objects that were attached above - it doesn't go off and do
    # the calculation again		
    foreach my $con (@{$new_vf->get_all_TranscriptVariations}) {
      
      my %gene_res;
      
      foreach my $string (@{$con->consequence_type}) {
	
	next if ( $string eq "INTERGENIC");

	if($con->cdna_start && $con->cdna_end && $con->cdna_start > $con->cdna_end) {
	  ($con->{'cdna_start'}, $con->{'cdna_end'}) = ($con->{'cdna_end'}, $con->{'cdna_start'});
	}
	
	if($con->translation_start &&  $con->translation_end && $con->translation_start > $con->translation_end) {
	  ($con->{'translation_start'}, $con->{'translation_end'}) = ($con->{'translation_end'}, $con->{'translation_start'});
	}

	if ( $con->transcript ) {

	  my $gene = $ga->fetch_by_transcript_stable_id($con->transcript->stable_id);
	  
	  $gene_res{ external_name } = $gene->external_name;
	  
	  $gene_res{ stable_id} = $gene->stable_id;
	  $gene_res{ transcript_id} = $con->transcript->stable_id;

	  my $xref = $con->transcript->get_all_DBEntries('RefSeq_dna' );

	  
	  if ( $$xref[0] ) {
	  
	    $gene_res{ xref } = $$xref[0]->display_id;
	    $gene_res{ transcript_id} = join("/",$gene_res{ xref },$gene_res{ transcript_id});
	  }

	  $gene_res{ cpos } = "";
	  $gene_res{ ppos } = "";


	  $gene_res{ position } = $string;
	  $gene_res{ cpos } = "c.".$con->cdna_start if ( $con->cdna_start);

	  if ( $con->translation_start) {
	    my ( $old, $new ) = split("\/", $con->pep_allele_string);
	    
	    $new = $old if ( $string eq "SYNONYMOUS_CODING");

	    $old = one2three( $old );
	    $new = one2three( $new );
	    

	    $gene_res{ ppos } = "p.$old".$con->translation_start . " $new";
	  }

	  $gene_res{ rs_number } = $existing_vf || "";
	  push @res, \%gene_res;
	} 
      }
    }
  }

  print Dumper( \@res );
  
  return \@res;
}




# 
# 
# 
# Kim Brugger (04 Jun 2010)
sub print_oneliner {
  my ( $pos, $res ) = @_;
  
  my @line;
  push @line, "$pos", "$$res{ref_base}>$$res{alt_base}";
  push @line, $$res{base_dist}{score};
  push @line, $$res{base_dist}{total};      

      map { push @line, $$res{base_dist}{$_} if ($$res{base_dist}{$_})} ( 'A', 'C', 'G', 'T', 'N');
      push @line, $$res{callers};

      if ($$res{snp_effect} ) {
      
	foreach my $snp_effect ( @{$$res{snp_effect}} ) {
	  my @effect_line;
	  push @effect_line, "$$snp_effect{ external_name }/$$snp_effect{ stable_id }", "$$snp_effect{ transcript_id }";
	  push @effect_line, $$snp_effect{ position };
	  push @effect_line, $$snp_effect{ cpos };
	  push @effect_line, $$snp_effect{ ppos };
	  push @effect_line, $$snp_effect{ rs_number } if ( $$snp_effect{ rs_number });

	  print join("\t", @line, @effect_line) . "\n";
	}
      }
      else {
	print join("\t", @line) . "\n";
      }



}


# 
# 
# 
# Kim Brugger (02 Jun 2010)
sub one2three {
  my ( $aminoacid) = @_;
  
  my %trans = (A => 'Ala',
	       R => 'Arg',
	       N => 'Asn',
	       D => 'Asp',
	       C => 'Cys',
	       E => 'Glu',
	       Q => 'Gln',
	       G => 'Gly',
	       H => 'His',
	       I => 'Ile',
	       L => 'Leu',
	       K => 'Lys',
	       M => 'Met',
	       F => 'Phe',
	       P => 'Pro',
	       S => 'Ser',
	       T => 'Thr',
	       W => 'Trp',
	       Y => 'Tyr',
	       V => 'Val');

  return $trans{ $aminoacid } if ($trans{ $aminoacid });
  return $aminoacid;
}


# 
# 
# 
# Kim Brugger (28 Apr 2010)
sub readin_bed {
  my ($file) = @_;

  my %indels = ();
  open (my $in, $file) || die "Could not open '$file': $!\n";
  while(<$in>) {

    next if (/^\#/);
    chomp;
    
    my ($chr, $start, $end, $change) = split("\t");

    $change =~ /(.)(\w+):(\d+)\/(\d+)/;
    my ( $type, $variation, $support, $depth) = ( $1, $2, $3, $4);
    
    $indels{ $chr }{ $start } = { end       => $end,
				  type      => $type,
				  variation => $variation,
				  support   => $support,
				  change    => $change };
    
  }



  return \%indels;
}
