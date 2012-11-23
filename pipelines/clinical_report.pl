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
use lib '/software/lib/e62/ensembl/modules/';
use lib '/software/lib/e62/ensembl-compara/modules/';
use lib '/software/lib/bioperl/';

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 0;
  if ( $DYNAMIC_LIB_PATHS ) {
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules" && ! grep /^$path\/modules/, @INC);
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    push @INC, '/home/kb468/easih-toolbox/modules/';
  }
}


use EASIH;
use EASIH::SNPs;
use EASIH::Git;

use EASIH::ChiSquare;


use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
#use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my @argv = @ARGV;

my %opts;
#getopts('cb:B:hi:l:o:O:Q:s:S:T', \%opts);
getopts('cb:B:h:i:l:o:O:Q:v:s:S:T', \%opts);    #svvd 31-10-2011: si -> v (snps,indels -> variants)
#perldoc() if ( $opts{h});
usage() if ( $opts{h});

my $species     = $opts{S} || "human";
my $buffer_size = 1;
my $host        = 'mgpc17.medschl.cam.ac.uk';
my $user        = "easih_ro";

# For Sticklers: CJP added 'NM_033150', 'NM_080630', 'NM_080679'; HM suggested to add: 'NM_001844', 'NM_001854'
# For TS: CJP added 'NM_001162427', 'NM_001114382'
my @ref_seqs = ();

# Howard wanted to change the low depth value to 30:
my $depth_cutoff = 30;

#$host = "ensembldb.ensembl.org";
#$user = "anonymous"; 

EASIH::SNPs->New();

if ( $opts{ Q }  ) {

  $opts{ Q } =~ s/\.bam//;
  $opts{ Q } =~ s/\.vcf//;
  $opts{ Q } =~ s/\.snps.vcf//;
  $opts{ Q } =~ s/\.indels.vcf//;

  $opts{s} = "$opts{Q}.snps.vcf"   if ( -e "$opts{Q}.snps.vcf");
  $opts{i} = "$opts{Q}.indels.vcf" if ( -e "$opts{Q}.indels.vcf");
  $opts{v} = "$opts{Q}.vcf" if ( -e "$opts{Q}.vcf");
  $opts{b} = "$opts{Q}.bam" if ( -e "$opts{Q}.bam" );
  $opts{o} = "$opts{Q}.var_full.csv";
  $opts{O} = "$opts{Q}.var.csv";
}


#print Dumper( \%opts );
  
#usage() if ( !$opts{s}  && !$opts{i});
usage() if ( !$opts{v} && !$opts{s} && !$opts{i});

#my $exit_count = 10;
my $exit_count = 9;

my $snp_vcf      = $opts{s};
my $indel_vcf    = $opts{i};
my $snp_indel_vcf= $opts{v};
my $bam          = $opts{b};
my $basecount    = $opts{c} || 0;
my $from_36      = $opts{T} || 0;

my $baits        = $opts{B} || "";
my $leeway       = 30;

my $bait_regions = readin_bed( $baits, $leeway ) if ( $baits );
my $out          = $opts{o} || undef;
my $filtered_out = $opts{O} || undef;

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


my $DACstring = "";
my $ANstring  = "";
my $SMstring  = "";




############################

open (*STDOUT, "> $out") || die "Could not open '$out': $!\n" if ( $out );
open( my $filtered_out_fh, "> $filtered_out") || die "Could not write to '$filtered_out': $!\n" if ( $filtered_out );
  
# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);

# get variation adaptors
my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
my $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');
my $sa  = $reg->get_adaptor($species, 'core', 'slice');
my $ga  = $reg->get_adaptor($species, 'core', 'gene');
my $ta  = $reg->get_adaptor($species, 'core', 'transcript');

# Getting the MethodLinkSpeciesSet adaptor: 
my $mlssa = $reg->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet');
#Get constrained element method_list_species_set
my $ce_mlss =  $mlssa->fetch_by_method_link_type_species_set_name("GERP_CONSTRAINED_ELEMENT", "mammals");
#Get constrained_element adaptor
my $ce_adaptor = $reg->get_adaptor('Multi', 'compara', 'ConstrainedElement');


my %slice_hash = ();
my ($sth_dbsnp, $sth_pop);

my $coverage = exon_coverage($bam, $bait_regions);


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
readin_vcf( $snp_indel_vcf ) if ( $snp_indel_vcf);


my %grch37_remapping = ();

my ($gatk_only, $pileup_only, $both_agree, $both_disagree) = (0, 0, 0, 0);
my $counter = 10;

print_var_header();

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
#    if (length($res{$position}{ref_base}) == 1 && length ($res{$position}{alt_base}) == 1) {
    $res{$position}{AAF}       = AAF( $chr, $pos, $res{$position}{ref_base}, $res{$position}{alt_base});
#    }

    
    
    my $Echr = $chr;
    $Echr =~ s/chr//;
    my $Epos = $pos;
    ($Echr, $Epos) = remap($Echr, $pos, $pos) if ( $from_36 );
    
#	$res{$position}{grcH37} = "$Echr:$Epos";
    
    next if ( ! $Echr );
    $grch37_remapping{"$chr:$pos"} = "$Echr:$Epos";
    my $slice = fetch_slice($Echr);
    next if ( ! $slice);
    my $allele_string = "$res{$position}{ref_base}/$res{$position}{alt_base}";
    $allele_string =~ s/,.*//;

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
    for my $l (@$line) {
      $l ='' if not defined $l;
    }
    $return_string .= join("\t", @$line) . "\n"; 
  }

  return $return_string;
}



# 
# 
# 
# Kim Brugger (08 Jul 2010)
sub print_var_header {

  my @res;

  $0 =~ s/.*\///;
  my @header;
  push @header, [ "#EASIH Variation Report v1.20.1"];
  push @header, [ "#commandline: $0 @argv"];
  push @header, [ "#dbases: ". EASIH::SNPs::db_info()];
  push @header, [ "#script version: $version"];
  
  push @header, [ "#bait filtering with a leeway of: $leeway and $baits as the bait file"] if ($baits );
  
  push @header, ["$coverage"] if ( $coverage );
  
    
  push @res, @header;
  
  my @annotations = ('Position', 'Change', 'Score', 'Depth', 'Genotype', 'AAF', 'FP');
  
  push @annotations, ('gene', 'transcript', 'exon', 'Effect', 'Nucleotide pos', 'AA change');
  push @annotations, ('dbsnp');
  push @annotations, ('dbsnp version');
  push @annotations, ('HGMD');
  
    
  push @res, [@annotations];
    


  print $filtered_out_fh text_table( \@res ) if ( $filtered_out );
  print STDOUT text_table( \@res );
  
  
}


# Kim Brugger (08 Jul 2010)
sub print_results {
  my ( $mapping ) = @_;

#  print STDERR "MAPPING DUMP :: " . Dumper( $mapping );
  my %done;
  my (@res, @filtered_res);

  foreach my $name ( sort keys %$mapping ) {

#    print STDERR "$name\n";



#    print STDERR  Dumper( $$mapping{$name}{res} );

    my $first = 1;

    foreach my $effect ( @{$$mapping{$name}{res}} ) {

#      print STDERR "$done{ $$effect{'transcript_id'}}\n";

      next if ($done{ $$effect{'transcript_id'}});
      $done{ $$effect{'transcript_id'}}++;

      my $transcript_id = $$effect{ transcript_id };
      $transcript_id =~ s/\.\d+//;
      
      my @effect_line;

      push @effect_line, "$name";

      if ( $$effect{ strand } == 1 ) {
	push @effect_line, "$$mapping{$name}{ref_base}>$$mapping{$name}{alt_base}";
      }
      else {
	push @effect_line, revDNA($$mapping{$name}{ref_base}).">".revDNA($$mapping{$name}{alt_base});
      }

      if ( $$effect{ HGVSc } ) {
	
	#ENST00000544455.1:c.3215T>A

 	$$effect{ HGVSc } =~ s/.*://;
 	$$effect{ HGVSc } =~ s/\D\>.*//;
 	$$effect{ HGVSc } =~ s/del.*//;
 	$$effect{ HGVSc } =~ s/ins.*//;

	$$effect{ cpos } = $$effect{ HGVSc };

      }
      else {
	next;
      }


      push @effect_line, $$mapping{$name}{qual};
      push @effect_line, $$mapping{$name}{depth};      
      push @effect_line, $$mapping{$name}{genotype};
      push @effect_line, $$mapping{$name}{AAF};
      push @effect_line, $$mapping{$name}{FP};
      
      push @effect_line, $$effect{ gene_id } || "";
      
      push @effect_line, $$effect{ transcript_id } || "";
      push @effect_line, $$effect{ exon }     || "";
      push @effect_line, $$effect{ effect }   || "";
      push @effect_line, $$effect{ cpos }     || "";
      push @effect_line, $$effect{ ppos }     || "";

      
      push @effect_line, $$effect{ rs_number   } || "";
      
      push @effect_line, $$effect{ HGMD   } || "";
      push @res, [@effect_line];
      
      push @filtered_res, [ @effect_line] if ( $first-- > 0);
    }
  }

  print $filtered_out_fh text_table( \@filtered_res ) if ( $filtered_out );
  print STDOUT text_table( \@res );
#  print STDOUT text_table( \@res );
}




# 
# 
# 
# Kim Brugger (11 May 2011)
sub sort_effects {
  my ($in_effects) = @_;


  my (@ensembl, @non_ensembl, @part_ensembl);

  foreach my $effect (@$in_effects ) {

    $$effect{ effect } ||= "INTERGENIC";
    $$effect{ gene_id } ||= "";
    $$effect{ transcript_id } ||= "";

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

  @ensembl      = sort { $effects{ $$b{effect} } || 0 <=> $effects{ $$a{effect} } || 0} @ensembl;
  @part_ensembl = sort { $effects{ $$b{effect} } || 0 <=> $effects{ $$a{effect} } || 0} @part_ensembl;
  @non_ensembl  = sort { $effects{ $$b{effect} } || 0 <=> $effects{ $$a{effect} } || 0} @non_ensembl;

  
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

    if ( 1 ) {
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
 	$HGMD        = $result->{hgmd};
      }
    }
    
    if (! $vf->slice->sub_Slice($vf->start, $vf->end, $vf->strand) ) {
      print STDERR "should this be hg18? cannot get slice: ". $vf->start .":". $vf->end ."-". $vf->strand . "\n";
      next;
    }

    # get consequences
    # results are stored attached to reference VF objects
    # so no need to capture return value here

    my $all_trans_var = $vf->get_all_TranscriptVariations;

    #print "LENGTH: [", scalar @$all_trans_var, "]\n";
    if (scalar @$all_trans_var == 0) {
      push @{$$vars{ $name }{res}}, {};
    }
    
    foreach my $tv (@{$all_trans_var}) {
      
      if($tv->cds_start && $tv->cds_end && $tv->cds_start > $tv->cds_end) {
	($tv->{'cds_start'}, $tv->{'cds_end'}) = ($tv->{'cds_end'}, $tv->{'cds_start'});
      }
      
      if($tv->translation_start &&  $tv->translation_end && $tv->translation_start > $tv->translation_end) {
	($tv->{'translation_start'}, $tv->{'translation_end'}) = ($tv->{'translation_end'}, $tv->{'translation_start'});
      }

      foreach my $tva (@{$tv->get_all_alternate_TranscriptVariationAlleles}) {



	my %gene_res;
    
#    $gene_res{ effect } = "INTERGENIC";
	$gene_res{ name } = $name;
#	print STDERR "TID:: : " . $tv->transcript->stable_id . "\n";
	my $gene = ($tv->transcript ? $ga->fetch_by_transcript_stable_id($tv->transcript->stable_id) : undef);
	my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};


	my @e = map { $_->display_term } @{$tva->get_all_OverlapConsequences};

#	print STDERR Dumper (@e );
	
	@e = sort { $effects{ $b } <=> $effects{ $a }} @e;

	$gene_res{strand}   = $tv->transcript->strand();
	
	$gene_res{ effect } = $e[0];
	$gene_res{ hgvs_coding } = $tva->hgvs_coding if ($tva->hgvs_coding);

	# HGVS
	if (defined($tva->hgvs_coding)) {
	  $gene_res{ 'HGVSc' } = $tva->hgvs_coding;
	  $gene_res{ 'HGVSc' } =~ s/c\.-/c.1-/;
	  if ($gene_res{ 'HGVSc' } =~ /c\.\*/) {
	    my $cds_end = $tv->transcript->length;
	    $gene_res{ 'HGVSc' } =~ s/c\.\*/c.$cds_end+/;
	  }

	    

	}

	$gene_res{ transcript_id }  = $tv->transcript->stable_id;
	$gene_res{ gene_id       }  = $gene->stable_id;

 	$gene_res{ rs_number }   = $existing_vf || "";

	if(scalar @entries) {
	  my @effects = @{$tv->consequence_type};
	  $gene_res{ effect } =  $effects[0];
	  $gene_res{ HGMD } = "Y" if (grep /'HGMD_MUTATION'/, @e);
	  
	  $gene_res{ gene_id } = $entries[0]->display_id;
	
	  my $xrefs = $tva->transcript->get_all_DBEntries('RefSeq_dna' );

	  # ensure that this is one of the transcripts we are looking for...
	  if ( $$xrefs[0] ) {
	    foreach my $xref ( @$xrefs ) {
	      
	      my $xref_id = $xref->display_id;
	      $xref_id =~ s/\.\d+//;
	      my @grep_res = grep { /^$xref_id/ } @ref_seqs;

	      next if ( @grep_res == 0 );
#	      print STDERR   " $xref_id ==  '@ref_seqs' [[ '@grep_res' ]] $gene_res{ transcript_id }\n";
	      $gene_res{ transcript_id } = $xref->display_id;

	      $gene_res{exon}     = exon_number($tv);
	    }
	  }

 	  my $gene = $ga->fetch_by_transcript_stable_id($tva->transcript->stable_id);
	  
 	  $gene_res{ cpos } = "?";
 	  $gene_res{ ppos } = "";

          if ( $tv->cds_start ) { 
# 	    $gene_res{ cpos } = "c.".$tv->cds_start;
          }
          else {
#            $gene_res{ cpos } = "c.".exon_offset($tv);
	  }
	
	
	  if ( $tv->translation_start) {
	    my ( $old, $new ) = ("","");
	    if ($tv->pep_allele_string) {
	      ( $old, $new ) = split("\/", $tv->pep_allele_string);
	      
	      $new = $old if ( !$new || $new eq "");
	      $old = one2three( $old );
	      $new = one2three( $new );
	      
	      if ($gene_res{ effect } eq 'FRAMESHIFT_CODING') {
		$new = "FS";
	      }

	      $gene_res{ ppos } = "p.$old".$tv->translation_start . " $new";
	    }
	  }
	}


#	print STDERR Dumper ( \%gene_res );
      
	push @{$$vars{ $name }{res}}, \%gene_res if ( $gene_res{ exon });
      }
    }
  }
}


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
# Kim Brugger (17 May 2011)
sub revDNA {
  my ($dna) = @_;

  $dna =~ tr/[ACGT]/[TGCA]/;
  $dna = reverse($dna);
  return $dna;
}


# 
# 
# 
# Kim Brugger (22 May 2012)
sub exon_offset {
  my ($tv) = @_;

    # work out which exon or intron this variant falls in

    # ensure the keys exist so even if we don't fall in an exon 
    # or intron we'll only call this method once


  my $vf = $tv->variation_feature;    
  
  my $vf_start = $vf->start;
  my $vf_end   = $vf->end;
  
  my $strand = $tv->transcript->strand;
  
  my $tran = $tv->transcript;

  # CJP - this seems to get alternative transcripts each time, from speaking with Howard, this should always be the same one,
  # so BRCA1 should always be transcript BRCA1-001 ENST00000357654 I think from comparing his exon numbering with those in EnsEMBL.
  # Maybe we need to do the same with other genes like BRCA2, Sticklers, etc.

  # TODO: One or two EnsEMBL exon starts/stops are slightly different than the clinical exon starts/stops
  # so this may need to be taken into account.
  # Maybe this also happens with BRCA2 exons?

#  if ($tran->external_name =~ /^BRCA1/) {
#    $tran = $ta->fetch_by_stable_id('ENST00000357654');
#  }

  my $exons = $tran->{_variation_effect_feature_cache}->{exons} ||= $tran->get_all_Exons;
#  @$exons = reverse( @$exons);
  
  my $tot_exons = scalar(@$exons);

  my $exon_count = 0;
  my $exon_sum = $$exons[0]->end - $$exons[0]->start+1;

  my $c_offset=0;

  

  for(my $i=1;$i< $tot_exons;$i++) {

    my $prev_exon_offset = -1;
    my $next_exon_offset = -1;

    my $intron_start = $strand == 1 ? $$exons[$i-1]->end + 1 : $$exons[$i]->end + 1;
    my $intron_end   = $strand == 1 ? $$exons[$i]->start - 1 : $$exons[$i-1]->start - 1;
      
    if ( $vf_start >= $intron_start && $vf_end <= $intron_end ) {

      my $pre_intron_offset = abs($vf_start - $intron_start);
      my $post_intron_offset = abs($vf_end - $intron_end);

      if ( $pre_intron_offset < $post_intron_offset) { 
#	return $$exons[$i-1]->cdna_end . " + " . abs($vf_start - $intron_start);
	return "$exon_sum + " . abs($vf_start - $intron_start);
      }
      else {
#	return $$exons[$i]->cdna_start . " - " . abs($vf_end - $intron_end);
	$exon_sum += $$exons[$i]->end - $$exons[$i]->start+1;
	return "$exon_sum - " . abs($vf_end - $intron_end);
      }
    }
    $exon_sum += $$exons[$i]->end - $$exons[$i]->start+1;
    
  }
  
}



# 
# 
# 
# Kim Brugger (22 May 2012)
sub exon_number {
  my ($tv) = @_;

  my $vf = $tv->variation_feature;    
  
  my $vf_start = $vf->start;
  my $vf_end   = $vf->end;
  my $vf_chr = $vf->slice->seq_region_name();
# = $tv->transcript->stable_id();

  my $xref = $tv->transcript->get_all_DBEntries('RefSeq_dna' );
  my $transcript = $$xref[0]->display_id if ( $$xref[0] );

  return undef if ( ! $transcript );
  $transcript =~ s/\.\d+//;

#  print "$transcript\n";
#  exit;

  my $exon_nr = undef;
  foreach my $region ( @{$$bait_regions{ $vf_chr }} ) {
    my ($reg_start, $reg_end, $exon, $reg_id) = @$region;


#    print STDERR "next if ( $reg_id ne $transcript );\n";

#    next if ( $reg_id ne $transcript );

    if ( $vf_start >= $reg_start && $vf_start <= $reg_end && $vf_end >= $reg_start && $vf_end <= $reg_end ) {
      return $exon;
    }

    $exon_nr = "INTRON";

#    if ( $vf_start >= $reg_start - $leeway && $vf_start < $reg_start) { 
#      return"INTRON";
#    }

#    if ( $vf_end > $reg_end && $vf_end <= $reg_end NM_001099857+ $leeway) { 
#      return "INTRON";
#    }

  }

  return $exon_nr;
}




# 
# 
# 
# Kim Brugger (22 May 2012)
sub exon_number_ensembl {
  my ($tv) = @_;

    # work out which exon or intron this variant falls in

    # ensure the keys exist so even if we don't fall in an exon 
    # or intron we'll only call this method once


  my $vf = $tv->variation_feature;    
  
  my $vf_start = $vf->start;
  my $vf_end   = $vf->end;
  
  my $strand = $tv->transcript->strand;
  my $tran = $tv->transcript;

  # CJP - this seems to get alternative transcripts each time, from speaking with Howard, this should always be the same one,
  # so BRCA1 should always be transcript BRCA1-001 ENST00000357654 I think from comparing his exon numbering with those in EnsEMBL.
  # Maybe we need to do the same with other genes like BRCA2, Sticklers, etc.

  if ($tran->external_name =~ /^BRCA1/) {
    $tran = $ta->fetch_by_stable_id('ENST00000357654');
  }

  my $exons = $tran->{_variation_effect_feature_cache}->{exons} ||= $tran->get_all_Exons;

  my $tot_exons = scalar(@$exons);

  my $exon_count = 0;

  my $prev_exon;
  my $exon_nr = "";
  my $intron_nr = "";

  for my $exon (@$exons) {

    $exon_count++;

    if ($exon_count == 4 and $tran->external_name =~ /^BRCA1\-001/) {
      $exon_count++; # To deal with clinical exon numbering system for BRCA1
    }
    
    if (overlap($vf_start, $vf_end, $exon->start, $exon->end)) {
      $exon_nr = $exon_count;
      last;
    }
    
    if ($prev_exon) {
      my $intron_start = $strand == 1 ? $prev_exon->end + 1 : $exon->end + 1;
      my $intron_end   = $strand == 1 ? $exon->start - 1 : $prev_exon->start - 1;
      
      if ($prev_exon && overlap($vf_start, $vf_end, $intron_start, $intron_end)) {
	$intron_nr = sprintf "%d/%d", $exon_count - 1, $tot_exons - 1;
	last;
      }
    }
    
    $prev_exon = $exon;
  }
  return $exon_nr;
}



sub overlap {
  my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
   
  return ( ($f1_end >= $f2_start) and ($f1_start <= $f2_end) );
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
sub AAF_two_alleles {
  my ( $chr, $SNP_pos, $ref_base, $alt_base) = @_;

#  print STDERR " $chr, $SNP_pos, $ref_base, $alt_base\n";

  if ( ! $bam ) {
#    print STDERR "need a bam file for finding base distribution\n";
    return;
  }

#  my %base_stats = ( A => 0, C => 0, G => 0, T => 0, N => 0);
  my %base_stats = ($alt_base =>0, $ref_base => 0);
  my %qual_stats;
  my $total = 0;

  open (my $st_pipe, "$samtools view -F0x404 $bam $chr:$SNP_pos-$SNP_pos | ") || die "Could not open samtools pipe: $!";

  while(<$st_pipe>) {
    chomp;
    my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");

    ($sequence, $quality) = patch_alignment($sequence, $quality, $cigar);

    my @seq   = split("", $sequence );
    my @qual  = split("", $quality  );

    my $base_pos = $SNP_pos - $pos;
    my $base = $seq[  $base_pos ];

    my $bases = "";

    use List::Util qw[max];    

    
    for(my $i = 0; $i < max(length( $ref_base), length( $alt_base)); $i++) {
      $bases .= $seq[  $base_pos + $i ] || "";
    }

    my $qual = $qual[ $base_pos ];

    $bases =~ s/-//g;

#    if ( $SNP_pos == 32891451 ) {
#      print "$sequence ($ref/$alt)\n";
#    }

    if ( !$base ) {
      print "$SNP_pos $_\n \n";
      print STDERR "FAILED !!!! \n";
      exit;
    }

#    next if ( (ord( $qual ) - 33) < 17);

    $base_stats{ $bases }++;


    $total++;
  }

  die Dumper( \%base_stats );


#  if ($base_stats{$alt_base} == 0 || $base_stats{$ref_base} == 0) {
#    print "$base_stats{$alt_base} == 0 || $base_stats{$ref_base} == 0\n";
#  }


#  print STDERR "$chr:$SNP_pos\n";

  my $noise     = 0;
  
  
#  $noise     = evaluate_SNP($base_stats{$ref_base}, $base_stats{$alt_base}) if (length($ref_base) == 1 && length($alt_base) == 1);

  print "$chr:$SNP_pos ref: $base_stats{$ref_base}, alt: $base_stats{$alt_base}, $noise\n";


#  exit     if ( $SNP_pos == 32891451 );

  return (1, $noise) if ($base_stats{$ref_base} == 0);
  return (0, $noise) if ($base_stats{$alt_base} == 0);
  return (sprintf("%.04f", $base_stats{$alt_base}/($base_stats{$alt_base}+$base_stats{$ref_base})), $noise);
}



# 
# Now using a python program to calculate this
# 
# Kim Brugger (15 Nov 2012)
sub AAF {
  my ( $chr, $SNP_pos, $ref_base, $alt_base) = @_;

#  return "0" if ($SNP_pos != 153791228);

  my @alt_alleles = split(/,/, $alt_base);

  use List::Util qw[max];    
  my $max_allele_lenght = max(length($ref_base), map { length($_) } @alt_alleles);

  map { $_ = $_."-"x(length($ref_base) - length($_)) } @alt_alleles;
  

  my $bam_allele_freq = "/home/kb468/easih-toolbox/scripts/bam_allele_freq.py";

  my $freqs_line = `$bam_allele_freq $bam $chr $SNP_pos`;
  chomp $freqs_line;
  
  my ( $pos, @freq_entries ) = split("\t", $freqs_line);
  my %freqs;
  my $depth = 0;
  foreach my $freq ( @freq_entries ) {
    my ($allele, $count) = split(":", $freq);
    $freqs{ $allele } = $count;
    $depth += $count;
  }

#  print STDERR "@alt_alleles\t$freqs_line\t$depth\n";
#  print STDERR Dumper( \%freqs );

  return join(",", map { sprintf("%.04f", $freqs{$_}/$depth)} @alt_alleles); 
}


# 
# 
# 
# Kim Brugger (29 Apr 2010)
sub AAF_old {
  my ( $chr, $SNP_pos, $ref_base, $alt_base) = @_;

  print STDERR " $chr, $SNP_pos, $ref_base, $alt_base\n";

  my @alt_alleles = split(/,/, $alt_base);
  my @alt_alleles_lengths  = map { length($_) } @alt_alleles;

  if ( ! $bam ) {
    print STDERR "need a bam file for finding base distribution\n";
    die;
  }

  use List::Util qw[max];    
  my $max_allele_lenght = max(length($ref_base), map { length($_) } @alt_alleles);

  my %base_stats = ($ref_base => 0);
  map { $base_stats{ $_ } = 0 } @alt_alleles;

  my %qual_stats;
  my $total = 0;

  open (my $st_pipe, "$samtools view -F0x404 $bam $chr:$SNP_pos-$SNP_pos | ") || die "Could not open samtools pipe: $!";

  while(<$st_pipe>) {
    chomp;
    my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");

    ($sequence, $quality) = patch_alignment($sequence, $quality, $cigar);

    my @seq   = split("", $sequence );
    my @qual  = split("", $quality  );

    my $base_pos = $SNP_pos - $pos;
    my $base = $seq[  $base_pos ];

    my $bases = "";
    
    for(my $i = 0; $i < $max_allele_lenght; $i++) {
      $bases .= $seq[  $base_pos + $i ] || "";
    }

    my $qual = $qual[ $base_pos ];

    $bases =~ s/-//g;


    if ( !$base ) {
      print "$SNP_pos $_\n \n";
      print STDERR "FAILED !!!! \n";
      exit;
    }

#    next if ( (ord( $qual ) - 33) < 17);

    $base_stats{ $bases }++;

    $total++;
  }


 # print STDERR Dumper( \%base_stats );


  if ( @alt_alleles == 1 ) {
    return 1 if ($base_stats{$ref_base} == 0);
    return 0 if ($base_stats{$alt_base} == 0);
  }
  my $depth = $base_stats{$ref_base};
  map { $depth += $base_stats{$_} } @alt_alleles;


  print STDERR Dumper( \%base_stats );

  print STDERR  "$chr, $SNP_pos " .  $base_stats{$_} ."/" . $base_stats{$_} ." '$depth'  @alt_alleles \n"; 


  return join(",", map { sprintf("%.04f", $base_stats{$_}/($base_stats{$_}+$depth))} @alt_alleles); 



  return sprintf("%.04f", $base_stats{$alt_base}/($base_stats{$alt_base}+$base_stats{$ref_base}));
}



# 
# 
# 
# Kim Brugger (13 Jun 2012)
sub evaluate_SNP {
  my ( $ref_count, $alt_count ) = @_;

  my $l_ref_count = 10*log($ref_count + 1);
  my $l_alt_count = 10*log($alt_count + 1);

  my  ($chisquare, $degrees_of_freedom, $chip) =  EASIH::ChiSquare::chisquare_raw( [$l_ref_count, $l_alt_count] );
  
  return 1 if ( $chip*100 < 5 && $ref_count > $alt_count);
#  return 0 if ( $chip*100 < 5);
  return 0;
  
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

#  print STDERR "$genotype, $reference\n";

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
      
      $SNPs{$chr}{$pos}{genotype} = "HOMO" if ($info_hash{AN} == 1);
      $SNPs{$chr}{$pos}{genotype} = "HET" if ($info_hash{AN} == 2);
      $SNPs{$chr}{$pos}{genotype} = "UNKNOWN" if ($info_hash{AN} != 1 && $info_hash{AN} != 2);
    }
   elsif ( ($info_hash{AC} && $info_hash{AC}  =~ /^(\d+),\d+\z/ )) {
     my $indels = $1;
     
     $SNPs{$chr}{$pos}{genotype} = "HOMO" if ( $indels*100/$depth > 75 );
     $SNPs{$chr}{$pos}{genotype} = "HET" if ( $indels*100/$depth <= 75 && $indels*100/$depth > 35 );
     $SNPs{$chr}{$pos}{genotype} = "UNKNOWN" if ( $indels*100/$depth <= 35 );
   }
  }


}




# 
# 
# 
# Kim Brugger (11 May 2010)
sub readin_bed {
  my ( $infile, $leeway ) = @_;

  my %res;

  my %ref_seqs;

  open (STDIN, $infile) || die "Could not open '$infile': $!\n" if ( $infile );
  while(<STDIN>) {

    chomp;
    next if (/^\z/);

    next if ((/^chromosome/) && (! $leeway)); 
    
    chomp;
    my ($chr, $start, $end, $name, $ref_id) = split("\t", $_);
    $ref_seqs{ $ref_id }++;

    ($chr, $start, $end) = $_ =~ /(.*?):(\d+)-(\d+)/
	if ( ! $start );

    next if ( ! $chr );

    $chr =~ s/chr//;
    
    
    push @{$res{$chr}}, [$start, $end, $name, $ref_id] if ( $chr);
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

  @ref_seqs = keys %ref_seqs;
  return \%res;
}
  



# 
# 
# 
# Kim Brugger (13 Jul 2010)
sub find_program {
  my ($program) = @_;


  my @paths = ("/home/easih/bin/",
	       "/home/kb468/bin/",
	       "/home/kb468/easih-toolbox/scripts/",
	       "/usr/local/bin");
  
  foreach my $path ( @paths ) {
    
    return "$path/$program" if ( -e "$path/$program" );
  }

  my $location = `which $program`;
  chomp( $location);
  
  return $location if ( $location );

  return undef;
}


sub DepthAndCoverage {
  my($bamfile, $bait_regions) = @_;

  my $samtools     = find_program('samtools');
  my $bam2depth    = find_program('bam2depth');

  my ($START, $END) = (0, 1);
  
  my $total_reads;
  my %base_coverage;
  my %exon_coverage;
  
  my %res;
  
  my $patched_start = 0;
  
  foreach my $chr ( sort keys %$bait_regions ) {
    
#  $chr ="chr2";
	
    my @regions =  @{$$bait_regions{$chr}};
   
    foreach my $region ( @regions ) {
      my (undef, $base_pos, $depth);
      my ($start, $end, $name) = @{$region};
 
      $name = "$chr:$start-$end" if ( ! $name );

      open (my $pipe, "$samtools depth -r $chr:$start-$end $bamfile  | " ) || die "Could not open pipe: $!\n";
      my ($min, $max, $sum) = (0,0,0);
      my $covered = 0;
      my $uncaptured = 0;
      while(<$pipe>) {
      
	
	(undef, $base_pos, $depth) = split("\t", $_);
	
	$min = $depth if ( $depth < $min && $min == 0);
	$max = $depth if ( $depth > $max );
	$covered++ if ( $depth > $depth_cutoff);
	
      }
      close( $pipe );
      
      if ( $covered < $end - $start + 1 ) {
	$uncaptured = $end - $start + 1 - $covered;
      }

      $res{ $name }{'min'} = $min;
      $res{ $name }{'max'} = $max;
      $res{ $name }{'uncaptured'} = $uncaptured;
      
      

    }
    
#  last;
  }

#  print Dumper( \%res );
  
      
#  return($DACstring);
}


sub exon_coverage {
  my($bamfile, $bait_regions) = @_;

  my $samtools     = find_program('samtools');
  my $bam2depth    = find_program('bam2depth');

  my ($START, $END) = (0, 1);
  
  
  my %res;
  
  my $patched_start = 0;
  
  foreach my $chr ( sort keys %$bait_regions ) {
    
#  $chr ="chr2";
	
    my @regions =  @{$$bait_regions{$chr}};

    foreach my $region ( @regions ) {
      my (undef, $base_pos, $depth);
      my ($start, $end, $name, $ref_id) = @{$region};
 
      $start -= $leeway;
      $end   += $leeway;

      $name = "$chr:$start-$end" if ( ! $name );

      my @depths;
      open (my $pipe, "$samtools depth -r $chr:$start-$end $bamfile  | " ) || die "Could not open pipe: $!\n";
      while(<$pipe>) {
	chomp;
	(undef, $base_pos, $depth) = split("\t", $_);
	$depths[ $base_pos - $start ] = $depth;
      }

      my ($min, $max, $sum) = (-1,0,0);
      my $covered = 0;
      my $uncaptured = 0;
      for( my $i = 0; $i < $end-$start+1;$i++) {
	
	my $depth = $depths[$i];
	$depth ||= 0;
	
	push @{$res{ $name }{'uncaptured'}}, $start + $i if ($depth == 0);
	push @{$res{ $name }{'low'}}, $start + $i if ($depth <= $depth_cutoff);

	$sum += $depth;
	$min = $depth if ( $depth < $min || $min == -1);
	$max = $depth if ( $depth > $max );
	$covered++;
      }
	
      close( $pipe );



      my @lows    = numbers2ranges($chr, $res{ $name }{'low'}, $ref_id );
      my @missing = numbers2ranges($chr, $res{ $name }{'uncaptured'}, $ref_id );

      
      $res{ $name }{'min'}  = $min;
      $res{ $name }{'max'}  = $max;
      $res{ $name }{'mean'} = sprintf("%.2f",$sum/($end-$start+1));
      $res{ $name }{'low'}  = \@lows;
      $res{ $name }{'missing'} = \@missing;
#      print Dumper( $res{ 'BRCA2_exon24' } ) if ( $name eq "BRCA2_exon24");

    }
    
#  last;
  }

  my @res_strings;
  push @res_strings, "#         ". join("\t", 'name', 'min', 'max', 'mean', 'lows', 'missing');

  
  foreach my $name ( sort {my ($a_pre, $a_post) = split("_exon", $a);
			   my ($b_pre, $b_post) = split("_exon", $b);
			   $a_pre cmp $b_pre || $a_post <=> $b_post
		     } keys %res ) {
    push @res_strings, "#  capture: ". join("\t", $name, $res{$name}{'min'}, $res{$name}{'max'}, $res{$name}{'mean'}, 
					   "'".join(',', @{$res{ $name }{'low'}})."'", 
					   "'".join(',', @{$res{ $name }{'missing'}})."'");
    
  }

  return join("\n", @res_strings);

#  return($DACstring);
}



# 
# 
# 
# Kim Brugger (29 May 2012)
sub numbers2ranges {
  my ($chr, $numbers, $ref_id ) = @_;
  
  my @res;
  my ($start, $end) = (-1,-1);
  foreach my $number (@$numbers) {
    if ( $start == -1 ) {
      $start = $number;
      $end   = $number;
      next;
    }
    elsif ( $number > $end + 1) {
      if ( $start == $end ) {
	push @res, coords2cpos($chr, $start, $ref_id);
      }
      else {
	push @res, coords2cpos($chr, $start, $ref_id)."-".coords2cpos($chr, $end, $ref_id);
      }
      
      ($start, $end) = ( $number, $number);
      
    }
    else {
      $end   = $number;
    }
  }	
  
  if ( $start != -1 ) {
    if ( $start == $end ) {
      push @res, coords2cpos($chr, $start, $ref_id);
    }
    else {
      push @res, coords2cpos($chr, $start, $ref_id)."-".coords2cpos($chr, $end, $ref_id);
    }
  }
  

#  print Dumper (\@res );
#  exit;
  return @res;
}



# 
# 
# 
# Kim Brugger (21 Aug 2012)
sub coords2cpos {
  my ($chr, $pos, $ref_id) = @_;

#  return $pos;

#  print "$chr, $pos, $ref_id c2p\n";

  my $slice = fetch_slice($chr);


  my %results;

  my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
    -start          => $pos,
    -end            => $pos,
    -slice          => $slice,           # the variation must be attached to a slice
    -allele_string  => "A/T",
    -strand         => 1,
    -map_weight     => 1,
    -adaptor        => $vfa,           # we must attach a variation feature adaptor
    -variation_name => $pos, # original position is used as the key!
      );

  my $effects = variation_effects(\%results, [$new_vf]);
 
  foreach my $effect ( @{$results{ $pos }{res}} ) {

#    print Dumper( $effect );

    $$effect{transcript_id} =~ s/\.\d+//;

    next if ($$effect{transcript_id} ne $ref_id);

    if ( $$effect{ HGVSc } ) {
      	
      $$effect{ HGVSc } =~ s/.*://;
      $$effect{ HGVSc } =~ s/\D\>.*//;
      $$effect{ HGVSc } =~ s/del.*//;
      $$effect{ HGVSc } =~ s/ins.*//;
      
      return $$effect{ HGVSc };
    }

  }
    

  return $pos;
}




# 
# 
# 
# Kim Brugger (09 Nov 2010)
sub usage {
  
  $0 =~ s/.*\///;

#  print "USAGE: $0 -b[am file] -i[indel vcf file] -s[np vcf file] -T<ranform, use if mapped against hg18> -B[ait file] -l[eeway, default 100 bp] -c[ount bases, need a -b as well]\n";
  print "USAGE: $0 -b[am file] -v[ariant vcf file] -T<ranform, use if mapped against hg18> -B[ait file] -l[eeway, default 100 bp] -c[ount bases, need a -b as well]\n";

#  print "\nor extrapolate the standard <bam, SNP vcf, indel vcf, output files> with the -Q <basefile name> option\n";
  print "\nor extrapolate the standard <bam, vcf, output files> with the -Q <basefile name> option\n";
  print "EXAMPLE: $0 -Q [base name] -T<ransform>\n";
  print "\n";

  
  print "USAGE: -o[output file]\n";  

  exit;
}



