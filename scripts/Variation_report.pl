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

use lib '/software/lib/e62/ensembl-variation/modules/';
use lib '/software/lib/e62/ensembl-functgenomics/modules/';
use lib '/software/lib/e62/ensembl/modules/';
use lib '/software/lib/e62/ensembl-compara/modules/';
use lib '/software/lib/bioperl/';

use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::SNPs;
use EASIH::Git;

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
#use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my @argv = @ARGV;

my %opts;
getopts('b:B:hi:l:o:O:Q:s:S:T', \%opts);
perldoc() if ( $opts{h});

my $species     = $opts{S} || "human";
my $buffer_size = 1;
my $host        = 'mgpc17';
my $user        = "easih_ro";

#$host = "ensembldb.ensembl.org";
#$user = "anonymous"; 

EASIH::SNPs->New();

if ( $opts{ Q }  ) {

  $opts{ Q } =~ s/\.bam//;
  $opts{ Q } =~ s/\.snps.vcf//;
  $opts{ Q } =~ s/\.indels.vcf//;

  $opts{s} = "$opts{Q}.snps.vcf";
  $opts{i} = "$opts{Q}.indels.vcf";
#  $opts{b} = "$opts{Q}.bam" if ( -e "$opts{Q}.bam" );
  $opts{o} = "$opts{Q}.var_full.csv";
  $opts{O} = "$opts{Q}.var.csv";
}
  

my $exit_count = 10;

my $snp_vcf      = $opts{s};
my $indel_vcf    = $opts{i};
my $bam          = $opts{b};
my $from_36      = $opts{T} || 0;
#my $at_ENSEMBL   = $opts{E} || 0;

my $baits        = $opts{B} || "";
my $leeway       = $opts{l} || 100;

my $bait_regions = readin_bed( $baits, $leeway ) if ( $baits );
my $out          = $opts{o} || undef;
my $full_out     = $opts{O} || undef;

my %effects  = ('ESSENTIAL_SPLICE_SITE'  => 10, 
		'STOP_GAINED'            => 10, 
		'STOP_LOST'              => 10, 
		'COMPLEX_INDEL'          => 10,
		'FRAMESHIFT_CODING'      => 10, 
		'NON_SYNONYMOUS_CODING'  => 10, 
		'HGMD_MUTATION'          =>  0, 
		'SPLICE_SITE'            =>  8, 
		'PARTIAL_CODON'          =>  8,
		'SYNONYMOUS_CODING'      =>  5, 
		'REGULATORY_REGION'      =>  3, 
		'WITHIN_MATURE_miRNA'    =>  3, 
		'5PRIME_UTR'             =>  1,
		'3PRIME_UTR'             =>  1, 
		'UTR'                    =>  1, 
		'INTRONIC'               =>  1, 
		'NMD_TRANSCRIPT'         =>  1, 
		'WITHIN_NON_CODING_GENE' =>  1, 
		'UPSTREAM'               =>  2,
		'DOWNSTREAM'             =>  2, 
		'INTERGENIC'             =>  2,
		'NO_CONSEQUENCE'         =>  0, 
    );

open (*STDOUT, "> $out") || die "Could not open '$out': $!\n" if ( $out );
open( my $full_out_fh, "> $full_out") || die "Could not write to '$full_out': $!\n" if ( $full_out );
  
# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);

# get variation adaptors
my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
my $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');
my $sa  = $reg->get_adaptor($species, 'core', 'slice');
my $ga  = $reg->get_adaptor($species, 'core', 'gene');

# Getting the MethodLinkSpeciesSet adaptor: 
my $mlssa = $reg->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet');
#Get constrained element method_list_species_set
my $ce_mlss =  $mlssa->fetch_by_method_link_type_species_set_name("GERP_CONSTRAINED_ELEMENT", "mammals");
#Get constrained_element adaptor
my $ce_adaptor = $reg->get_adaptor('Multi', 'compara', 'ConstrainedElement');


my %slice_hash = ();
my ($sth_dbsnp, $sth_pop);
my $use_local_dbsnp = 1;#$opts{L} || 0;

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

my $version = EASIH::Git::version();


my %SNPs = ();

my $samtools = `which samtools`;
chomp( $samtools);

#readin_pileup( $pileup ) if ( $pileup);
readin_vcf( $snp_vcf ) if ( $snp_vcf);
readin_vcf( $indel_vcf ) if ( $indel_vcf);

my %grch37_remapping = ();

my ($gatk_only, $pileup_only, $both_agree, $both_disagree) = (0, 0, 0, 0);
my $counter = 10;
foreach my $chr ( sort {$a cmp $b}  keys %SNPs ) {
  
  my %res;
  my @vfs;

  foreach my $pos ( sort { $a <=> $b} keys %{$SNPs{$chr}} ) {
        
    my $position = "$chr:$pos";
    
    $res{$position}{ref_base}  = $SNPs{$chr}{$pos}{ref_base};
    $res{$position}{alt_base}  = $SNPs{$chr}{$pos}{alt_base};
    $res{$position}{qual}      = $SNPs{$chr}{$pos}{qual};
    $res{$position}{depth}     = $SNPs{$chr}{$pos}{depth};    
    $res{$position}{filter}    = $SNPs{$chr}{$pos}{filter};
    $res{$position}{genotype}  = $SNPs{$chr}{$pos}{genotype};
    $res{$position}{base_dist} = base_dist( $chr, $pos, $res{$position}{ref_base}, $res{$position}{alt_base}) if ( $bam );
    
    my $Echr = $chr;
    $Echr =~ s/chr//;
    my $Epos = $pos;
    ($Echr, $Epos) = remap($Echr, $pos, $pos) if ( $from_36 );
    
#	$res{$position}{grcH37} = "$Echr:$Epos";
    
    next if ( ! $Echr );
    $grch37_remapping{"$chr:$pos"} = "$Echr:$Epos";
    my $slice = fetch_slice($Echr);
      my $allele_string = "$res{$position}{ref_base}/$res{$position}{alt_base}";
      
      # create a new VariationFeature object
      my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
	-start          => $Epos,
	-end            => $Epos,
	-slice          => $slice,           # the variation must be attached to a slice
	-allele_string  => $allele_string,
	-strand         => 1,
	-map_weight     => 1,
	-adaptor        => $vfa,           # we must attach a variation feature adaptor
	-variation_name => $position, # original position is used as the key!
	  );
      
      push @vfs, $new_vf; 
      
      
      if ( @vfs >= $buffer_size) {
	my $effects = variation_effects(\%res, \@vfs);
	print_results( \%res, $effects);
	undef @vfs;# = ();	  
	undef %res;# = ();
	%grch37_remapping = ();
#	  exit if ($exit_count-- < 0);
      }
    }
  
  
  
  if ( @vfs || keys %res ) {
    my $effects = variation_effects(\@vfs);
    print_results( \%res, $effects);
    @vfs = ();
    %res = ();
    %grch37_remapping = ();
  }

}



#
# Creates a simple table, the function expects an array of arrays, and a border or not flag.
# 
# Kim Brugger (20 Oct 2003)
sub text_table {
  my ($cells) = @_;


  my $return_string = "";
  foreach my $line ( @$cells ) {
    $return_string .= join("\t", @$line) . "\n";
  }

  return $return_string;
}


my ($printed_header) = (0);

# 
# 
# 
# Kim Brugger (08 Jul 2010)
sub print_results {
  my ( $mapping ) = @_;

#  print Dumper($effects);

  my (@res, @filtered_res);
  if ( ! $printed_header++ ) {
    $0 =~ s/.*\///;
    my @header;
    push @header, [ "#EASIH Variation Report v1.20"];
    push @header, [ "#commandline: $0 @argv"];
    push @header, [ "#dbases: ". EASIH::SNPs::db_info()];
    push @header, [ "#version: $version"];

    push @header, [ "#bait filtering with a leeway of: $leeway and $baits as the bait file"] if ($baits );

    push @res, @header;
    push @filtered_res, @header;
    push @filtered_res, ["#filtered for most important snp effect"];
    
    my @annotations = ('Position', 'Change', 'Filter', 'Score', 'Depth', 'Genotype');
    push @annotations, ('', '', '', '', '') if ( $bam );
      
    push @annotations, ('gene', 'transcript', 'Effect', 'codon pos', 'AA change');
    push @annotations, ('Grantham score');
    push @annotations, ('dbsnp');
    push @annotations, ('dbsnp flags');
    push @annotations, ('HGMD');

    push @annotations, 'pfam';
    push @annotations, 'PolyPhen';
    push @annotations, 'SIFT';
    push @annotations, 'Condel';
    push @annotations, 'GERP';

    push @res, [@annotations];
    push @filtered_res, [@annotations];
  }

  foreach my $name ( sort keys %$mapping ) {

    my @line;
    push @line, "$name";
    
    push @line, "$$mapping{$name}{ref_base}>$$mapping{$name}{alt_base}";
    push @line, $$mapping{$name}{filter};
    push @line, $$mapping{$name}{qual};
    push @line, $$mapping{$name}{depth};      
    push @line, $$mapping{$name}{genotype};
    if ( $$mapping{$name}{base_dist} ) {
      map { push @line, $$mapping{$name}{base_dist}{$_} if ($$mapping{$name}{base_dist}{$_})} ( 'A', 'C', 'G', 'T', 'N');
    }


#    @{$$mapping{$name}{res}} = sort { $effects{ $$b{effect} } <=> $effects{ $$a{effect} } } @{$$mapping{$name}{res}};

    $$mapping{$name}{res} = sort_effects( $$mapping{$name}{res});

    my $first = 1;

    foreach my $effect ( @{$$mapping{$name}{res}} ) {

    # find the most important/interesting variation effect
#    foreach my $snp_effect ( @$effect ) {
#    my $snp_effect = $$mapping{ $name }{ res };
      my @effect_line;
	
      push @effect_line, $$effect{ gene_id } || "";
      
      push @effect_line, $$effect{ transcript_id } || "";
      push @effect_line, $$effect{ effect }   || "";
      push @effect_line, $$effect{ cpos }     || "";
      push @effect_line, $$effect{ ppos }     || "";
      push @effect_line, $$effect{ grantham } || "";
      
      push @effect_line, $$effect{ rs_number   } || "";
      push @effect_line, $$effect{ dbsnp_flags } || "" if ( $use_local_dbsnp );
      
      push @effect_line, $$effect{ HGMD   } || "";
      push @effect_line, $$effect{ pfam } || "";
      
      foreach my $tool (qw(SIFT PolyPhen Condel)) {
	my $lc_tool = lc($tool);
	push @effect_line, $$effect{$tool} || "";
      }

      push @effect_line, $$effect{ gerp } || "";
      
      push @res, [@line, @effect_line];
      
      push @filtered_res, [@line, @effect_line] if ( $first-- > 0);
    }
  }

  print $full_out_fh text_table( \@filtered_res          ) if ( $full_out );
  print STDOUT text_table( \@res );
}




# 
# 
# 
# Kim Brugger (11 May 2011)
sub sort_effects {
  my ($in_effects) = @_;


  my (@ensembl, @non_ensembl, @part_ensembl);

  foreach my $effect (@$in_effects ) {

    $$effect{ gene_id} ||= "INTERGENIC";

    if ( $$effect{ gene_id}       !~ /^ENSG\d+/ &&
	 $$effect{ transcript_id} !~ /^ENST\d+/ ) {
      push @non_ensembl, $effect;
    }
    elsif ( $$effect{ gene_id}       !~ /^ENSG\d+/ &&
	    $$effect{ transcript_id} =~ /^ENST\d+/ ) {
      push @part_ensembl, $effect;
    }
    else {
      push @ensembl, $effect;
    }

  }    

  @ensembl      = sort { $effects{ $$b{effect} } <=> $effects{ $$a{effect} }} @ensembl;
  @part_ensembl = sort { $effects{ $$b{effect} } <=> $effects{ $$a{effect} }} @part_ensembl;
  @non_ensembl  = sort { $effects{ $$b{effect} } <=> $effects{ $$a{effect} }} @non_ensembl;

  
  return [@non_ensembl, @part_ensembl, @ensembl];
}



# 
# 
# 
# Kim Brugger (26 Apr 2011), contact: kim.brugger@easih.ac.uk
sub variation_effects {
  my ($vars, $var_features) = @_;

  foreach my $vf (@$var_features) {    

    my $name = $vf->variation_name();

    my $existing_vf = "";
    my ($dbsnp_flags, $HGMD) = ('','', '');
    if ( $use_local_dbsnp ) {
      my ($chr, $pos);
      if ($from_36 && $grch37_remapping{$vf->variation_name()}) {
	($chr, $pos) = split(":", $grch37_remapping{$vf->variation_name()});
	#	print "looking for snp at: $chr $pos ($existing_vf) ". ($vf->variation_name())."\n";
      }
      else {
 	($chr, $pos) = split(":", $vf->variation_name());
      }
      
      my $result = EASIH::SNPs::fetch_snp($chr, $pos);
      
      if ( $result->{rs} ) {
 	$existing_vf = $result->{rs};
 	$dbsnp_flags = $result->{flags} || "";
	$dbsnp_flags .= ";VC=$result->{class}" if ($result->{class});
 	$HGMD        = $result->{hgmd};
      }
    }


    my $cons = $ce_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($ce_mlss,$vf->slice->sub_Slice($vf->start, $vf->end, $vf->strand));

    # get consequences
    # results are stored attached to reference VF objects
    # so no need to capture return value here
    foreach my $tv (@{$vf->get_all_TranscriptVariations}) {

      if($tv->cdna_start && $tv->cdna_end && $tv->cdna_start > $tv->cdna_end) {
	($tv->{'cdna_start'}, $tv->{'cdna_end'}) = ($tv->{'cdna_end'}, $tv->{'cdna_start'});
      }
      
      if($tv->translation_start &&  $tv->translation_end && $tv->translation_start > $tv->translation_end) {
	($tv->{'translation_start'}, $tv->{'translation_end'}) = ($tv->{'translation_end'}, $tv->{'translation_start'});
      }

      foreach my $tva (@{$tv->get_all_alternate_TranscriptVariationAlleles}) {

	my %gene_res;
    
#    $gene_res{ effect } = "INTERGENIC";
	$gene_res{ name } = $name;	
	my $gene = ($tv->transcript ? $ga->fetch_by_transcript_stable_id($tv->transcript->stable_id) : undef);
	my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};

	my @e = map { $_->display_term } @{$tva->get_all_OverlapConsequences};
	@e = sort { $effects{ $b } <=> $effects{ $a }} @e;
	
	$gene_res{ effect } = $e[0];
	$gene_res{hgvs_coding} = $tva->hgvs_coding if ($tva->hgvs_coding);
	$gene_res{gerp} = $$cons[0]->score if ( @$cons );
	$gene_res{ transcript_id }  = $tv->transcript->stable_id;
	$gene_res{ gene_id       }  = $gene->stable_id;


 	$gene_res{ rs_number }   = $existing_vf || "";
 	$gene_res{ dbsnp_flags } = $dbsnp_flags || "";

	if(scalar @entries) {
	  my @effects = @{$tv->consequence_type};
	  $gene_res{ effect } =  $effects[0];
	  $gene_res{ HGMD } = "Y" if (grep /'HGMD_MUTATION'/, @e);
	  
	  $gene_res{ gene_id } = $entries[0]->display_id;
	
	  my $xref = $tva->transcript->get_all_DBEntries('RefSeq_dna' );
	  $gene_res{ transcript_id } = $$xref[0]->display_id if ( $$xref[0] );
 	  my $gene = $ga->fetch_by_transcript_stable_id($tva->transcript->stable_id);
	  
 	  $gene_res{ cpos } = "";
 	  $gene_res{ ppos } = "";

 	  $gene_res{ cpos } = "c.".$tv->cdna_start if ( $tv->cdna_start);

 	  if ( $tv->translation_start) {
 	    my ( $old, $new ) = ("","");
 	    if ($tv->pep_allele_string) {
 	      ( $old, $new ) = split("\/", $tv->pep_allele_string);
	      
 	      $new = $old if ( !$new || $new eq "");
 	      $old = one2three( $old );
 	      $new = one2three( $new );
 	      $gene_res{ ppos } = "p.$old".$tv->translation_start . " $new";
 	      $gene_res{ grantham } = grantham_score($old, $new);
 	    }

 	    my $protein = $tv->transcript->translation();

 	    my $prot_feats = $protein->get_all_ProteinFeatures();
	    
 	    while (my $prot_feat = shift @{ $prot_feats }) {
 	      my $logic_name = $prot_feat->analysis()->logic_name();
	      
 	      next if ( $logic_name ne 'Pfam');
	      
 	      if ($tva->translation_start >= $prot_feat->start() and
 		  $tva->translation_end <= $prot_feat->end() ) {
		
 		$gene_res{ pfam }     = $prot_feat->idesc();
 		$gene_res{ interpro } = $prot_feat->interpro_ac();
 	      } 
 	    }
 	  }
 	}
	  

	# HGVS
	$gene_res{ 'HGVSc' } = $tva->hgvs_coding if defined($tva->hgvs_coding);
	$gene_res{ 'HGVSp' } = $tva->hgvs_protein if defined($tva->hgvs_protein);
	
	foreach my $tool (qw(SIFT PolyPhen Condel)) {
	  my $lc_tool = lc($tool);
	  
	  my $pred_meth   = $lc_tool.'_prediction';
	  my $score_meth  = $lc_tool.'_score';
	  
	  my $pred = $tva->$pred_meth;
	  my $score = $tva->$score_meth;
	  
	  $gene_res{$tool} = "$pred" if (!$score && $pred);
	  $gene_res{$tool} = "$score" if ($score && !$pred);
	  $gene_res{$tool} = "$pred/$score" if (defined $score && defined $pred);
	}
	
	push @{$$vars{ $name }{res}}, \%gene_res;
      }
    }
  }
}


#     my $name = $vf->variation_name();
#     # find any co-located existing VFs
#     my $existing_vf = "";
#     my $regulations;
    
#     if(defined($vf->adaptor->db)) {
#       my $fs = $vf->feature_Slice;

#       if($fs->start > $fs->end) {
# 	($fs->{'start'}, $fs->{'end'}) = ($fs->{'end'}, $fs->{'start'});
#       }
#       foreach my $existing_vf_obj(@{$vf->adaptor->fetch_all_by_Slice($fs)}) {
# 	$existing_vf = $existing_vf_obj->variation_name
# 	    if ($existing_vf_obj->seq_region_start == $vf->seq_region_start &&
# 		$existing_vf_obj->seq_region_end   == $vf->seq_region_end );
#       }

     
#     }



# #	$existing_vf = "";
#     my ($dbsnp_flags, $dbsnp_freq, $HGMD, $phylop, $phast) = ('','', '', '', '');
#     if ( $use_local_dbsnp ) {
#       my ($chr, $pos);
#       if ($from_36 && $grch37_remapping{$vf->variation_name()}) {
# 	($chr, $pos) = split(":", $grch37_remapping{$vf->variation_name()});
# #	print "looking for snp at: $chr $pos ($existing_vf) ". ($vf->variation_name())."\n";
#       }
#       else {
# 	($chr, $pos) = split(":", $vf->variation_name());
#       }

#       my $result = EASIH::SNPs::fetch_snp_GRCh37($chr, $pos);

#       if ( $result->{rs} ) {
# 	$existing_vf = $result->{rs};
# 	$dbsnp_flags = $result->{flags} || "";
# 	$dbsnp_freq  = EASIH::SNPs::population_stats($existing_vf);
# 	$HGMD        = $result->{hgmd};
#       }

#       $phylop = EASIH::SNPs::phylop_score_GRCh37($chr, $pos); 
#       $phast  = EASIH::SNPs::phast_score_GRCh37($chr, $pos); 
	  
#     }

	
#     # the get_all_TranscriptVariations here now just retrieves the
#     # objects that were attached above - it doesn't go off and do
#     # the calculation again		
#     foreach my $tva (@{$vf->get_all_TranscriptVariations}) {

#       my %gene_res;
#       $gene_res{ name } = $name;
      
#       foreach my $string (@{$tva->consequence_type}) {
	
# 	$gene_res{ position } = $string;

# 	if($tva->cdna_start && $tva->cdna_end && $tva->cdna_start > $tva->cdna_end) {
# 	  ($tva->{'cdna_start'}, $tva->{'cdna_end'}) = ($tva->{'cdna_end'}, $tva->{'cdna_start'});
# 	}
	
# 	if($tva->translation_start &&  $tva->translation_end && $tva->translation_start > $tva->translation_end) {
# 	  ($tva->{'translation_start'}, $tva->{'translation_end'}) = ($tva->{'translation_end'}, $tva->{'translation_start'});
# 	}

# 	if ( $con->transcript ) {


# 	  my $gene = $ga->fetch_by_transcript_stable_id($con->transcript->stable_id);
	  
# 	  $gene_res{ external_name } = $gene->external_name;
# 	  $gene_res{ stable_id}      = $gene->stable_id;
# 	  $gene_res{ transcript_id}  = $con->transcript->stable_id;

# 	  my $xref = $con->transcript->get_all_DBEntries('RefSeq_dna' );
  
# 	  $gene_res{ xref } = $$xref[0]->display_id
# 	      if ( $$xref[0] );

# 	  $gene_res{ cpos } = "";
# 	  $gene_res{ ppos } = "";


# 	  $gene_res{ cpos } = "c.".$con->cdna_start if ( $con->cdna_start);



# 	  if ( $con->translation_start) {
# 	    my ( $old, $new ) = ("","");
# 	    if ($con->pep_allele_string) {
# 	      ( $old, $new ) = split("\/", $con->pep_allele_string);
	      
# 	      $new = $old if ( !$new || $new eq "");
# 	      $old = one2three( $old );
# 	      $new = one2three( $new );
# 	      $gene_res{ ppos } = "p.$old".$con->translation_start . " $new";
# 	      $gene_res{ grantham } = grantham_score($old, $new);
# 	    }

# 	    my $protein = $con->transcript->translation();

# 	    my $prot_feats = $protein->get_all_ProteinFeatures();
	    
# 	    while (my $prot_feat = shift @{ $prot_feats }) {
# 	      my $logic_name = $prot_feat->analysis()->logic_name();
	      
# 	      next if ( $logic_name ne 'Pfam');
	      
# 	      if ($con->translation_start >= $prot_feat->start() and
# 		  $con->translation_end <= $prot_feat->end() ) {
		
# 		$gene_res{ pfam }     = $prot_feat->idesc();
# 		$gene_res{ interpro } = $prot_feat->interpro_ac();
# 	      } 
# 	    }
# 	  }
# 	}

# 	$gene_res{ regulation } = $regulations; 

# 	$gene_res{ rs_number }   = $existing_vf || "";
# 	$gene_res{ dbsnp_freq }  = $dbsnp_freq || "";
# 	$gene_res{ dbsnp_flags } = $dbsnp_flags || "";
# 	$gene_res{ HGMD   } = $HGMD || "";
# 	$gene_res{ phylop } = $phylop || "";
# 	$gene_res{ phast  } = $phast || "";

# 	push @{$res[$feature]}, \%gene_res;

#       }
#     }
#     $feature++;
#   }

#   return \@res;
# }


sub format_coords {
	my ($start, $end) = @_;
	
	if(!defined($start)) {
		return '-';
	}
	elsif(!defined($end)) {
		return $start;
	}
	elsif($start == $end) {
		return $start;
	}
	elsif($start > $end) {
		return $end.'-'.$start;
	}
	else {
		return $start.'-'.$end;
	}
}


# 
# 
# 
# Kim Brugger (28 May 2010)
sub variation_effects_old {
  my ($var_features) = @_;


  my @res = ();
  my $feature = 0;
  
  # get consequences
  # results are stored attached to reference VF objects
  # so no need to capture return value here
  $tva->fetch_all_by_VariationFeatures( $var_features );
  foreach my $vf (@$var_features) {    

    my $name = $vf->variation_name();
    # find any co-located existing VFs
    my $existing_vf = "";
    
    if(defined($vf->adaptor->db)) {
      my $fs = $vf->feature_Slice;

      if($fs->start > $fs->end) {
	($fs->{'start'}, $fs->{'end'}) = ($fs->{'end'}, $fs->{'start'});
      }
      foreach my $existing_vf_obj(@{$vf->adaptor->fetch_all_by_Slice($fs)}) {
	$existing_vf = $existing_vf_obj->variation_name
	    if ($existing_vf_obj->seq_region_start == $vf->seq_region_start &&
		$existing_vf_obj->seq_region_end   == $vf->seq_region_end );
      }

      
    }



#	$existing_vf = "";
    my ($dbsnp_flags, $dbsnp_freq, $HGMD, $phylop, $phast) = ('','', '', '', '');
    if ( $use_local_dbsnp ) {
      my ($chr, $pos);
      if ($from_36 && $grch37_remapping{$vf->variation_name()}) {
	($chr, $pos) = split(":", $grch37_remapping{$vf->variation_name()});
#	print "looking for snp at: $chr $pos ($existing_vf) ". ($vf->variation_name())."\n";
      }
      else {
	($chr, $pos) = split(":", $vf->variation_name());
      }

      my $result = EASIH::SNPs::fetch_snp_GRCh37($chr, $pos);

      if ( $result->{rs} ) {
	$existing_vf = $result->{rs};
	$dbsnp_flags = $result->{flags} || "";
	$dbsnp_freq  = EASIH::SNPs::population_stats($existing_vf);
	$HGMD        = $result->{hgmd};
      }

      $phylop = EASIH::SNPs::phylop_score_GRCh37($chr, $pos); 
      $phast  = EASIH::SNPs::phast_score_GRCh37($chr, $pos); 
	  
    }

	
    # the get_all_TranscriptVariations here now just retrieves the
    # objects that were attached above - it doesn't go off and do
    # the calculation again		
    foreach my $con (@{$vf->get_all_TranscriptVariations}) {

      my %gene_res;
      $gene_res{ name } = $name;
      
      foreach my $string (@{$con->consequence_type}) {
	
	$gene_res{ position } = $string;

	if($con->cdna_start && $con->cdna_end && $con->cdna_start > $con->cdna_end) {
	  ($con->{'cdna_start'}, $con->{'cdna_end'}) = ($con->{'cdna_end'}, $con->{'cdna_start'});
	}
	
	if($con->translation_start &&  $con->translation_end && $con->translation_start > $con->translation_end) {
	  ($con->{'translation_start'}, $con->{'translation_end'}) = ($con->{'translation_end'}, $con->{'translation_start'});
	}

	if ( $con->transcript ) {


	  my $gene = $ga->fetch_by_transcript_stable_id($con->transcript->stable_id);
	  
	  $gene_res{ external_name } = $gene->external_name;
	  $gene_res{ stable_id}      = $gene->stable_id;
	  $gene_res{ transcript_id}  = $con->transcript->stable_id;

	  my $xref = $con->transcript->get_all_DBEntries('RefSeq_dna' );
  
	  $gene_res{ xref } = $$xref[0]->display_id
	      if ( $$xref[0] );

	  $gene_res{ cpos } = "";
	  $gene_res{ ppos } = "";


	  $gene_res{ cpos } = "c.".$con->cdna_start if ( $con->cdna_start);



	  if ( $con->translation_start) {
	    my ( $old, $new ) = ("","");
	    if ($con->pep_allele_string) {
	      ( $old, $new ) = split("\/", $con->pep_allele_string);
	      
	      $new = $old if ( !$new || $new eq "");
	      $old = one2three( $old );
	      $new = one2three( $new );
	      $gene_res{ ppos } = "p.$old".$con->translation_start . " $new";
	      $gene_res{ grantham } = grantham_score($old, $new);
	    }

	    my $protein = $con->transcript->translation();

	    my $prot_feats = $protein->get_all_ProteinFeatures();
	    
	    while (my $prot_feat = shift @{ $prot_feats }) {
	      my $logic_name = $prot_feat->analysis()->logic_name();
	      
	      next if ( $logic_name ne 'Pfam');
	      
	      if ($con->translation_start >= $prot_feat->start() and
		  $con->translation_end <= $prot_feat->end() ) {
		
		$gene_res{ pfam }     = $prot_feat->idesc();
		$gene_res{ interpro } = $prot_feat->interpro_ac();
	      } 
	    }
	  }
	}


	$gene_res{ rs_number }   = $existing_vf || "";
	$gene_res{ dbsnp_freq }  = $dbsnp_freq || "";
	$gene_res{ dbsnp_flags } = $dbsnp_flags || "";
	$gene_res{ HGMD   } = $HGMD || "";
	$gene_res{ phylop } = $phylop || "";
	$gene_res{ phast  } = $phast || "";

	push @{$res[$feature]}, \%gene_res;

      }
    }
    $feature++;
  }

  return \@res;
}



# 
# 
# 
# Kim Brugger (06 Jul 2010)
sub remap {
  my ($chr, $start, $end, $strand) = @_;
  
  $strand ||= 0;

  my @res = $mapper->map($chr, $start, $end, $strand, $cs_from);
  foreach my $res ( @res ) {
    if ( $res->isa( 'Bio::EnsEMBL::Mapper::Coordinate' )) {
      my $chr_slice = $sa->fetch_by_seq_region_id($res->id);
      return ($chr_slice->seq_region_name, $res->start, $res->end, $res->strand);
    }
  }

  print STDERR "Could not remap: $chr, $start, $end\n";

  return (undef, undef, undef, undef);
}



# 
# cached slice fetcher..
# 
# Kim Brugger (06 Jul 2010)
sub fetch_slice {
  my ( $chr) = @_;

  my $slice;
  # check if we have fetched this slice already
  # disabled the caching, as we have a memory problem on our hands...
  if( defined $slice_hash{$chr}) {
    $slice = $slice_hash{$chr};
  }
 
  # if not create a new one
  else {
    
    # first try to get a chromosome
    eval { $slice = $sa->fetch_by_region('chromosome', $chr); };

    eval { $slice = $sa->fetch_by_region('supercontig', $chr); } if(!defined($slice));
    
    # if failed, try to get any seq region
    if(!defined($slice)) {
      $slice = $sa->fetch_by_region(undef, $chr);
    }
    
    # if failed, die
    if(!defined($slice)) {
      die("ERROR: Could not fetch slice named $chr\n");
    }	
    # store the hash
    $slice_hash{$chr} = $slice;
  }

  return $slice;
}



# 
# 
# 
# Kim Brugger (02 Jun 2010)
sub one2three {
  my ( $aminoacid) = @_;

  return "" if ( ! defined $aminoacid);
  
  my %trans = ('A' => 'Ala',
	       'R' => 'Arg',
	       'N' => 'Asn',
	       'D' => 'Asp',
	       'C' => 'Cys',
	       'E' => 'Glu',
	       'Q' => 'Gln',
	       'G' => 'Gly',
	       'H' => 'His',
	       'I' => 'Ile',
	       'L' => 'Leu',
	       'K' => 'Lys',
	       'M' => 'Met',
	       'F' => 'Phe',
	       'P' => 'Pro',
	       'S' => 'Ser',
	       'T' => 'Thr',
	       'W' => 'Trp',
	       'Y' => 'Tyr',
	       'V' => 'Val',
	       '*' => 'Ter');

  return $trans{ $aminoacid } if ($trans{ $aminoacid });
  return $aminoacid;
}


# 
# 
# 
# Kim Brugger (29 Apr 2010)
sub base_dist {
  my ( $chr, $SNP_pos, $ref_base, $alt_base) = @_;

  if ( ! $bam ) {
#    print STDERR "need a bam file for finding base distribution\n";
    return;
  }

  my %base_stats = ( A => 0, C => 0, G => 0, T => 0, N => 0);
  my %qual_stats;
  my $total = 0;

  open (my $st_pipe, "$samtools view -F0x404 $bam $chr:$SNP_pos-$SNP_pos | ") || die "Could not open samtools pipe: $!";

  while(<$st_pipe>) {
    chomp;
    my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");

    ($sequence, $quality) = patch_alignment($sequence, $quality, $cigar);

    my @seq  = split("", $sequence);
    my @qual = split("", $quality);    

    my $base_pos = $SNP_pos - $pos;
    my $base = $seq[ $base_pos ];

    if ( !$base ) {
      print "$SNP_pos $_\n \n";
      print STDERR "FAILED !!!! \n";
      exit;
    }

    my $qual = ord($qual[ $base_pos ])-33;
    $base_stats{ $base }++;
    push @{$qual_stats{$base}},  $qual;
    $total++;
  }

  foreach my $key ( keys %qual_stats ) {
    
    my $sum = eval join '+', @{$qual_stats{$key}};
    my $count = @{$qual_stats{$key}};
    my $avg_qual      =  sprintf("%.2f",$sum / $count);
    $avg_qual ||= 0;
    $qual_stats{$key} = $avg_qual;
  }

  my %res;
  my $alt_perc = 0;
  foreach my $base (sort {$base_stats{$b} <=> $base_stats{$a}} keys %base_stats ) {
    next if ( ! $total );
    my $perc = sprintf("%.2f", $base_stats{$base}/$total*100);
    $alt_perc = $perc if ( $base eq $alt_base );
    my $qual = $qual_stats{$base} || 0;
    $res{$base} = "$base: $base_stats{$base}($perc%)/$qual";
  }

  $res{type}  = "homozygous"   if ( $alt_perc > 76);
  $res{type}  = "heterozygous" if ( $alt_perc < 75 && $alt_perc > 35);
  $res{type}  = "unknown"      if ( $alt_perc < 35);
  $res{type} .= " ($alt_perc)";

  $res{total} = $total;
  
  return \%res;
}


# 
# 
# 
# Kim Brugger (20 Jul 2009)
sub patch_alignment {
  my ( $seq, $qual, $cigar ) = @_;

  return ($seq, $qual) if ( $cigar !~ /[DIS]/);
  
  my @seq  = split("", $seq );
  my @qual = split("", $qual );


  my (@cigar) = $cigar =~ /(\d*\w)/g;

  my $offset = 0;

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)



  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d+)(\w)/;

    if ( $type eq 'M') {
      $offset += $length;
      next;
    }
    elsif ( $type eq "D") {
      my @dashes = split("", "-"x$length);
      splice(@seq,  $offset, 0, @dashes);
      splice(@qual, $offset, 0, @dashes);
      $offset += $length;
    }
    elsif ( $type eq "I" || $type eq "S" ) {
      splice(@seq,  $offset, $length);
      splice(@qual, $offset, $length);
    }    

  }


  return (join("", @seq), join("",@qual));
}



# 
# 
# 
# Kim Brugger (29 Apr 2010)
sub genotype_2_base {
  my ($genotype) = @_;

  my ($base1, $base2) = split("", $genotype);
  return $base1 if ( $base1 eq $base2);
  
  if ( $base1 gt $base2 ) {
    my $tbase = $base1;
    $base1 = $base2;
    $base2 = $tbase;
  }

  return 'W' if ( $base1 eq 'A' && $base2 eq 'T' );
  return 'S' if ( $base1 eq 'C' && $base2 eq 'G' );
  return 'M' if ( $base1 eq 'A' && $base2 eq 'C' );
  return 'K' if ( $base1 eq 'G' && $base2 eq 'T' );
  return 'R' if ( $base1 eq 'A' && $base2 eq 'G' );
  return 'Y' if ( $base1 eq 'C' && $base2 eq 'T' );

  print STDERR "Does not know how to translate '$genotype'\n";

  return $base1;
}


# 
# 
# 
# Kim Brugger (29 Apr 2010)
sub subtract_reference {
  my ($genotype, $reference) = @_;

  print STDERR "$genotype, $reference\n";

  return $genotype if ( $genotype eq 'A' ||
			$genotype eq 'C' ||
			$genotype eq 'G' ||
			$genotype eq 'T' );

  my ($base1, $base2);



  ($base1, $base2) = ('A', 'T') if ( $genotype eq 'W');
  ($base1, $base2) = ('C', 'G') if ( $genotype eq 'S');
  ($base1, $base2) = ('A', 'C') if ( $genotype eq 'M');
  ($base1, $base2) = ('G', 'T') if ( $genotype eq 'K');
  ($base1, $base2) = ('A', 'G') if ( $genotype eq 'R');
  ($base1, $base2) = ('C', 'T') if ( $genotype eq 'Y');

  return $base1 if ( $reference eq $base2);
  return $base2 if ( $reference eq $base1);
  return $genotype;
  
#  print"Cannot subtract $reference from $genotype ( $base1, $base2) \n";
#  exit -1;
}



# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub in_bait_region_simple {
  my ($chr, $pos, $baits) = @_;

  die "I should not be here!!! \n";

  $baits ||= $bait_regions;

  $chr =~ s/chr//;
  
  return 0 if ( $chr ne "X" && $pos != 57421888);
  
  foreach my $region ( @{$$baits{$chr}})  {
    my ( $start, $end) = @$region;

    print STDERR "$start >= $pos && $end <= $pos\n";

    return 1 if ( $start <= $pos && $end >= $pos);
    last if ( $end > $pos );
  }

  print STDERR "Discarding $chr $pos\n";


#  $dropped++;  
  return 0;
  
}


# 
# 
# 
# Kim Brugger (31 Mar 2010)
sub verbose {
  return;
  my ($message, $level) = @_;
  $level  = 1;
  $message =~ s/\n+\Z//g;
  print STDERR "MESS " . ":"x$level . " $message\n";
#  sleep 1;
}


# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub in_bait_region {
  my ($chr, $pos, $baits) = @_;

  $baits ||= $bait_regions;
  $chr =~ s/chr//;

  return if ( ! $$baits{$chr});

  use POSIX qw(ceil floor);


  my $START = 0;
  my $END   = 1;
  
#  return 0 if ( $pos != 57421888);

  my @regions = @{$$baits{$chr}};

  # set the start and end of the array and find the 
  # the middle of the array
  my ( $left, $right ) = (0, int(@regions));
  my $middle = floor(($right - $left)/2);

  # Flush the buffer constantly
  $| = 1;

  my $loop_counter = 0;
    
  while (1) {

    verbose("MIDDLE $middle ( $left, $right)\n", 1);
    verbose(" $pos <=> $regions[ $middle ][$START] $regions[ $middle ][$END]\n", 1);
    
    # The new block is to the left of the middle.
    if ( $pos < $regions[ $middle ][$START] ) {
      $right = $middle;
      $middle = $left + floor(($right - $left)/2);
      verbose("L");
      last if ( $right <= $left || $middle == $left || $middle == $right);
    }
    # The new block is to the right of the middle.
    elsif ($pos > $regions[ $middle ][$END] ) {
      $left = $middle;
      $middle = $left + floor(($right - $left)/2);
      verbose("R");
      last if ( $right <= $left || $middle == $left || $middle == $right);
    }
    #
    # Now things gets interesting, we here start to calculate
    # overlapping and contained regions.
    #
    
    # this is a contained snp, exactly what we want!!!!
    elsif ( $pos >= $regions[ $middle ][ $START ]  &&
	    $pos <= $regions[ $middle ][ $END   ] ) {
      
      verbose("CONTAINED BLOCK", 2);
      return 1;
      last;
    }
    else {
      last;
    }
  }

  verbose("Discarding $chr $pos\n");


#  $dropped++;  
  return 0;
}



# 
# 
# 
# Kim Brugger (28 Apr 2010)
sub readin_vcf {
  my ($file) = @_;
  open (my $in, $file) || die "Could not open '$file': $!";

  my $used = 0;
  my $dropped = 0;

  while(<$in>) {
    next if (/^\#/);
    
    my ($chr, $pos, $id, $ref_base, $alt_base, $qual, $filter, $info) = split("\t");

    if ($bait_regions && ! in_bait_region($chr, $pos)) {
      $dropped++;
      next;
    }

    $used++;

#    $alt_base = subtract_reference($alt_base, $ref_base);

    my %info_hash;
    foreach my $entry (split("\;", $info )) {
      my @f=split("\=", $entry); 
      $info_hash{$f[0]} = $f[1];
    }
    

    my $depth = $info_hash{ DP };

    
    $SNPs{$chr}{$pos} = { depth     => $depth,
			  qual      => $qual,
			  alt_base  => $alt_base,
			  filter    => $filter,
			  pos       => $pos,
			  ref_base  => $ref_base};


    # for snps it should be one number
    if ($info_hash{AF} ) {
      
      $SNPs{$chr}{$pos}{genotype} = "HOMO" if ($info_hash{AF} == 1);
      $SNPs{$chr}{$pos}{genotype} = "HET" if ($info_hash{AF} == 0.50);
      $SNPs{$chr}{$pos}{genotype} = "UNKNOWN" if ($info_hash{AF} != 1 && $info_hash{AF} != 0.5);
    }
    elsif ( ($info_hash{AC} && $info_hash{AC}  =~ /^(\d+),\d+\z/ )) {
      my $indels = $1;
      
      $SNPs{$chr}{$pos}{genotype} = "HOMO" if ( $indels*100/$depth > 75 );
      $SNPs{$chr}{$pos}{genotype} = "HET" if ( $indels*100/$depth <= 75 && $indels*100/$depth > 35 );
      $SNPs{$chr}{$pos}{genotype} = "UNKNOWN" if ( $indels*100/$depth <= 35 );
    }
    else {
      $SNPs{$chr}{$pos}{genotype} = "unknown $info_hash{IAC}";
    }


  }

  print STDERR "Used: $used, Dropped: $dropped\n" if ($bait_regions);

}



 
# 
# 
# 
# Kim Brugger (28 Apr 2010)
sub readin_pileup {
  my ($file, $min_SNP_qual) = @_;
  open (my $in, $file) || die "Could not open '$file': $!\n";

  die "Pileup is no more\n";

  $min_SNP_qual ||= 0;

  while(<$in>) {
    next if (/^\#/);

    my ($chr, $pos, $ref_base, $alt_base, $cons_qual, $cons_SNP_qual, $max_map_qual, $depth, $pile, $quals) = split("\t");

#    print "$_";
    $alt_base = subtract_reference($alt_base, $ref_base);

    next if ( $cons_SNP_qual < $min_SNP_qual);
    
    $SNPs{$chr}{$pos}{samtools} = { depth        => $depth,
				    mapping_qual => $cons_SNP_qual,
				    alt_base     => $alt_base,
				    pos => $pos,};

    $SNPs{$chr}{$pos}{ref_base} = $ref_base;
  }
  
}
 


# 
# 
# 
# Kim Brugger (11 May 2010)
sub readin_bed {
  my ( $infile, $leeway ) = @_;

  my %res;

  open (STDIN, $infile) || die "Could not open '$infile': $!\n" if ( $infile );
  while(<STDIN>) {

    chomp;
    my ($chr, $start, $end) = split("\t", $_);

    ($chr, $start, $end) = $_ =~ /(.*?):(\d+)-(\d+)/
	if ( ! $start );

    next if ( ! $chr );

    $chr =~ s/chr//;
    
    $start -= $leeway;
    $end   += $leeway;
    
    push @{$res{$chr}}, [$start, $end] if ( $chr);
  }

  foreach my $key ( keys %res ) {
    
    @{$res{$key}} = sort { $$a[0] <=> $$b[0] } @{$res{$key}};
    my @tmp;
    my @data = @{$res{$key}};
    
    for(my $i=0;$i< @data; $i++) {
      
      # need at least one element in the array, so push and move on.
      if ( ! @tmp ) {
	push @tmp, $data[ $i ];
	next;
      }
      
      # contained in the region
      if ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ]  &&
	   $data[ $i ][ 1 ] <= $tmp[ -1 ][ 1 ] ) {
	next;
      }
      # overlapping
      elsif ( $data[ $i ][ 0 ] >= $tmp[ -1 ][ 0 ]  &&
	      $data[ $i ][ 0 ] <= $tmp[ -1 ][ 1 ]) {
	
	$tmp[ -1 ][ 1 ] = $data[ $i ][ 1 ];
      }
      # There is a gap between the end block and this one. Just push it on the end of  the array!
      else {
	push @tmp, $data[ $i ];
      }
    }
    @{$res{$key}} = @tmp;
  }

  return \%res;
}
  


# 
# 
# 
# Kim Brugger (01 Dec 2010)
sub grantham_score {
  my ($aa1, $aa2) = @_;

  $aa1 = uc($aa1);
  $aa2 = uc($aa2);

  return 0 if ( $aa1 eq $aa2);
  
  ($aa1, $aa2) = ($aa2, $aa1) if ($aa1 gt $aa2);
  
  my %grantham;
  $grantham{ALA}{ARG}=112; 
  $grantham{ALA}{ASN}=111; 
  $grantham{ALA}{ASP}=126; 
  $grantham{ALA}{CYS}=195; 
  $grantham{ALA}{GLN}=91; 
  $grantham{ALA}{GLU}=107; 
  $grantham{ALA}{GLY}=60; 
  $grantham{ALA}{HIS}=86; 
  $grantham{ALA}{ILE}=94; 
  $grantham{ALA}{LEU}=96; 
  $grantham{ALA}{LYS}=106; 
  $grantham{ALA}{MET}=84; 
  $grantham{ALA}{PHE}=113; 
  $grantham{ALA}{PRO}=27; 
  $grantham{ALA}{SER}=99; 
  $grantham{ALA}{THR}=58; 
  $grantham{ALA}{TRP}=148; 
  $grantham{ALA}{TYR}=112; 
  $grantham{ALA}{VAL}=64;
  
  $grantham{ARG}{ASN}=86; 
  $grantham{ARG}{ASP}=96; 
  $grantham{ARG}{CYS}=180; 
  $grantham{ARG}{GLN}=43; 
  $grantham{ARG}{GLU}=54; 
  $grantham{ARG}{GLY}=125; 
  $grantham{ARG}{HIS}=29; 
  $grantham{ARG}{ILE}=97; 
  $grantham{ARG}{LEU}=102; 
  $grantham{ARG}{LYS}=26; 
  $grantham{ARG}{MET}=91; 
  $grantham{ARG}{PHE}=97; 
  $grantham{ARG}{PRO}=103; 
  $grantham{ARG}{SER}=110; 
  $grantham{ARG}{THR}=71; 
  $grantham{ARG}{TRP}=101; 
  $grantham{ARG}{TYR}=77; 
  $grantham{ARG}{VAL}=96; 

  $grantham{ASN}{ASP}=23; 
  $grantham{ASN}{CYS}=139; 
  $grantham{ASN}{GLN}=46; 
  $grantham{ASN}{GLU}=42; 
  $grantham{ASN}{GLY}=80; 
  $grantham{ASN}{HIS}=68; 
  $grantham{ASN}{ILE}=149; 
  $grantham{ASN}{LEU}=153; 
  $grantham{ASN}{LYS}=94; 
  $grantham{ASN}{MET}=142; 
  $grantham{ASN}{PHE}=158; 
  $grantham{ASN}{PRO}=91; 
  $grantham{ASN}{SER}=46; 
  $grantham{ASN}{THR}=65;
  $grantham{ASN}{TRP}=174;
  $grantham{ASN}{TYR}=143;
  $grantham{ASN}{VAL}=133;

  $grantham{ASP}{CYS}=154; 
  $grantham{ASP}{GLN}=61; 
  $grantham{ASP}{GLU}=45; 
  $grantham{ASP}{GLY}=94; 
  $grantham{ASP}{HIS}=81; 
  $grantham{ASP}{ILE}=168; 
  $grantham{ASP}{LEU}=172; 
  $grantham{ASP}{LYS}=101; 
  $grantham{ASP}{MET}=160; 
  $grantham{ASP}{PHE}=177; 
  $grantham{ASP}{PRO}=108; 
  $grantham{ASP}{SER}=65; 
  $grantham{ASP}{THR}=85; 
  $grantham{ASP}{TRP}=181; 
  $grantham{ASP}{TYR}=160; 
  $grantham{ASP}{VAL}=152;

  $grantham{CYS}{GLN}=154; 
  $grantham{CYS}{GLU}=170; 
  $grantham{CYS}{GLY}=159; 
  $grantham{CYS}{HIS}=174; 
  $grantham{CYS}{ILE}=198; 
  $grantham{CYS}{LEU}=198; 
  $grantham{CYS}{LYS}=202; 
  $grantham{CYS}{MET}=196; 
  $grantham{CYS}{PHE}=205; 
  $grantham{CYS}{PRO}=169; 
  $grantham{CYS}{SER}=112; 
  $grantham{CYS}{THR}=149; 
  $grantham{CYS}{TRP}=215; 
  $grantham{CYS}{TYR}=194; 
  $grantham{CYS}{VAL}=192;

  $grantham{GLN}{GLU}=29; 
  $grantham{GLN}{GLY}=87; 
  $grantham{GLN}{HIS}=24; 
  $grantham{GLN}{ILE}=109; 
  $grantham{GLN}{LEU}=113; 
  $grantham{GLN}{LYS}=53; 
  $grantham{GLN}{MET}=101; 
  $grantham{GLN}{PHE}=116; 
  $grantham{GLN}{PRO}=76;
  $grantham{GLN}{SER}=68; 
  $grantham{GLN}{THR}=42; 
  $grantham{GLN}{TRP}=130; 
  $grantham{GLN}{TYR}=99; 
  $grantham{GLN}{VAL}=96;

  $grantham{GLU}{GLY}=98; 
  $grantham{GLU}{HIS}=40; 
  $grantham{GLU}{ILE}=134; 
  $grantham{GLU}{LEU}=138; 
  $grantham{GLU}{LYS}=56; 
  $grantham{GLU}{MET}=126; 
  $grantham{GLU}{PHE}=140; 
  $grantham{GLU}{PRO}=93; 
  $grantham{GLU}{SER}=80; 
  $grantham{GLU}{THR}=65; 
  $grantham{GLU}{TRP}=152; 
  $grantham{GLU}{TYR}=122; 
  $grantham{GLU}{VAL}=121;

  $grantham{GLY}{HIS}=89; 
  $grantham{GLY}{ILE}=135; 
  $grantham{GLY}{LEU}=138; 
  $grantham{GLY}{LYS}=127; 
  $grantham{GLY}{MET}=127; 
  $grantham{GLY}{PHE}=153; 
  $grantham{GLY}{PRO}=42; 
  $grantham{GLY}{SER}=56; 
  $grantham{GLY}{THR}=59; 
  $grantham{GLY}{TRP}=184; 
  $grantham{GLY}{TYR}=147; 
  $grantham{GLY}{VAL}=109;
  
  $grantham{HIS}{ILE}=94; 
  $grantham{HIS}{LEU}=99; 
  $grantham{HIS}{LYS}=32; 
  $grantham{HIS}{MET}=87; 
  $grantham{HIS}{PHE}=100; 
  $grantham{HIS}{PRO}=77; 
  $grantham{HIS}{SER}=89; 
  $grantham{HIS}{THR}=47; 
  $grantham{HIS}{TRP}=115; 
  $grantham{HIS}{TYR}=83; 
  $grantham{HIS}{VAL}=84; 
  
  $grantham{ILE}{LEU}=5; 
  $grantham{ILE}{LYS}=102; 
  $grantham{ILE}{MET}=10; 
  $grantham{ILE}{PHE}=21; 
  $grantham{ILE}{PRO}=95; 
  $grantham{ILE}{SER}=142; 
  $grantham{ILE}{THR}=89; 
  $grantham{ILE}{TRP}=61; 
  $grantham{ILE}{TYR}=33; 
  $grantham{ILE}{VAL}=29;
  
  $grantham{LEU}{LYS}=107; 
  $grantham{LEU}{MET}=15; 
  $grantham{LEU}{PHE}=22; 
  $grantham{LEU}{PRO}=98; 
  $grantham{LEU}{SER}=145; 
  $grantham{LEU}{THR}=92; 
  $grantham{LEU}{TRP}=61; 
  $grantham{LEU}{TYR}=36; 
  $grantham{LEU}{VAL}=32;
  
  $grantham{LYS}{MET}=95; 
  $grantham{LYS}{PHE}=102; 
  $grantham{LYS}{PRO}=103; 
  $grantham{LYS}{SER}=121; 
  $grantham{LYS}{THR}=78; 
  $grantham{LYS}{TRP}=110; 
  $grantham{LYS}{TYR}=85; 
  $grantham{LYS}{VAL}=97;
  
  $grantham{MET}{PHE}=28; 
  $grantham{MET}{PRO}=87; 
  $grantham{MET}{SER}=135; 
  $grantham{MET}{THR}=81; 
  $grantham{MET}{TRP}=67; 
  $grantham{MET}{TYR}=36; 
  $grantham{MET}{VAL}=21;
  
  $grantham{PHE}{PRO}=114; 
  $grantham{PHE}{SER}=155; 
  $grantham{PHE}{THR}=103; 
  $grantham{PHE}{TRP}=40; 
  $grantham{PHE}{TYR}=22; 
  $grantham{PHE}{VAL}=50;
  
  $grantham{PRO}{SER}=74; 
  $grantham{PRO}{THR}=38; 
  $grantham{PRO}{TRP}=147; 
  $grantham{PRO}{TYR}=110; 
  $grantham{PRO}{VAL}=68;
  
  $grantham{SER}{THR}=58; 
  $grantham{SER}{TRP}=177; 
  $grantham{SER}{TYR}=144; 
  $grantham{SER}{VAL}=124;
  
  $grantham{THR}{TRP}=128; 
  $grantham{THR}{TYR}=92; 
  $grantham{THR}{VAL}=69;
  
  $grantham{TRP}{TYR}=37; 
  $grantham{TRP}{VAL}=88;
  
  $grantham{TYR}{VAL}=55;


  return $grantham{$aa1}{$aa2} if ($grantham{$aa1}{$aa2});
  return "NA";

}




# 
# 
# 
# Kim Brugger (09 Nov 2010)
sub usage {
  
  $0 =~ s/.*\///;

  print "USAGE: $0 -b[am file] -i[indel vcf file] -s[np vcf file] -T<ranform, use if mapped against hg18> -B[ait file] -l[eeway, default 10 bp]\n";

  print "\nor extrapolate the standard <bam, SNP vcf, indel vcf, output files> with the -Q <basefile name> option\n";
  print "EXAMPLE: $0 -Q [base name] -T<ransform>\n";
  print "\n";

  
  print "USAGE: -o[output file]\n";  
  print "USAGE: -O[output file, filtered, one line/snp]\n";  

  exit;
}




# 
# 
# 
# Kim Brugger (09 Jul 2010)
sub perldoc {

  system "perldoc $0";
  exit;
}


=pod

=head1 SYNOPSIS

SNP_report turns a vcf file into a nice report, annotated with Ensembl information like effect, regulation etc

=head1 OPTIONS


=over

=item B<-b F<bam file>>: 

The bamfile that the SNP calling was based on. This is used for doing the base distribution information for each snp.

=item B<-f>: 

Prints out a multi-line report, otherwise it is done on a oneline basis (good for excel)

=item B<-H>: 

Prints a HTML report, default is a tab-separated one.

=item B<-I>: 

Adds links to IGV to ease navigation.

=item B<-m>: 

Minumum mapping quality to use for the base distribution report


=item B<-p>: 

Extracts pfam domains information from ensembl

=item B<-q>: 

Minumum SNP quality to report

=item B<-v F<vcf file>>: 

The infile that is to be converted into a nice report.

=item B<-T>: 

The mapping was done on NCBI36/hg18 and the coordinates should be transformed to GRCh37/hg19.

=back

=head1 NOTES

=over

=item B<-s>: will change the species, but changing that will probably break about just everything.


=cut
