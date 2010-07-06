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

use lib '/home/kb468/projects/ensembl-variation/modules/';
use lib '/home/kb468/projects/e57/ensembl/modules/';
use lib '/home/kb468/projects/e57/bioperl-live/';

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

my %opts;
getopts('s:g:p:b:T', \%opts);

my $species     = $opts{s} || "human";
my $buffer_size = 1;
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
my %slice_hash = ();

my $gatk        = $opts{g};
my $pileup      = $opts{p};
my $bam         = $opts{b};
my $from_36     = $opts{T} || 0;
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


my %SNPs = ();

my $samtools = `which samtools`;
chomp( $samtools);

readin_pileup( $pileup ) if ( $pileup);
readin_cvf( $gatk )      if ( $gatk);

my %reports;

my ($gatk_only, $pileup_only, $both_agree, $both_disagree) = (0, 0, 0, 0);
my @vfs;

foreach my $chr ( sort {$a cmp $b}  keys %SNPs ) {

  # Hack so we only look at chrX for develop purposes
  $chr = 'chrX';


  my %res;

  foreach my $pos ( sort { $a <=> $b} keys %{$SNPs{$chr}} ) {
 
    
    my $position = "$chr:$pos";
    

    my @keys = keys %{$SNPs{$chr}{$pos}};

    @keys = grep(!/ref_base/, @keys);
    
    $res{$position}{callers} = join("/", @keys);

    if (1 || keys %{{ map {$SNPs{$chr}{$pos}{$_}{alt_base}, 1} @keys }} == 1) {
	
      my $key = "GATK";
      $res{$position}{ref_base} = $SNPs{$chr}{$pos}{ref_base};
      $res{$position}{alt_base} = $SNPs{$chr}{$pos}{$key}{alt_base};

      $res{$position}{base_dist}  = base_dist( $chr, $pos, $res{$position}{ref_base}, $res{$position}{alt_base});

      my $Echr = $chr;
      $Echr =~ s/chr//;
      ($Echr, $pos) = remap($Echr, $pos, $pos) if ( $from_36 );

      my $slice = fetch_slice($Echr);
      my $allele_string = "$res{$position}{ref_base}/$res{$position}{alt_base}";

      # create a new VariationFeature object
      my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
	-start          => $pos,
	-end            => $pos,
	-slice          => $slice,           # the variation must be attached to a slice
	-allele_string  => $allele_string,
	-strand         => 1,
	-map_weight     => 1,
	-adaptor        => $vfa,           # we must attach a variation feature adaptor
	-variation_name => $position, # original position is used as the key!
	  );

      push @vfs, $new_vf; 


    }
    else {

      print "BROKEN/NOT USED ANYMORE!!! '$res{$position}{callers}'\n";
      next;

      foreach my $key ( @keys ) {
	$res{ref_base} = $SNPs{$chr}{$pos}{ref_base};
	$res{alt_base} = $SNPs{$chr}{$pos}{$key}{alt_base};
	$SNPs{$chr}{$pos}{base_dist}  = base_dist( $chr, $pos, $res{ref_base}, $res{alt_base});
	$res{snp_effect} = snp_effect($chr, $pos, $pos, "$res{ref_base}/$res{alt_base}");
	print_oneliner( $position, \%res);
      }
    }

    if ($buffer_size <= @vfs ) {
      my $effects = variation_effects(\@vfs);
      print_oneliner( \%res, $effects);
      @vfs = ();
      %res = ();
    }

    if ( ! $full_report && ! $html_out ) {
    }      
#    exit;
    

  }

  if ( @vfs ) {
    my $effects = variation_effects(\@vfs);
    print_oneliner( \%res, $effects);
    @vfs = ();
    %res = ();
  }


  last;
}



# 
# 
# 
# Kim Brugger (04 Jun 2010)
sub print_oneliner {
  my ( $mapping, $effects ) = @_;

  foreach my $effect ( @$effects ) {

    my $name = $$effect[0]{ name };
    my @line;
    push @line, "$name", "$$mapping{$name}{ref_base}>$$mapping{$name}{alt_base}";
    if ( $$mapping{$name}{base_dist} ) {
      push @line, $$mapping{$name}{base_dist}{score};
      push @line, $$mapping{$name}{base_dist}{total};      
    
      map { push @line, $$mapping{$name}{base_dist}{$_} if ($$mapping{$name}{base_dist}{$_})} ( 'A', 'C', 'G', 'T', 'N');
      push @line, $$mapping{$name}{callers};
    }
    
    foreach my $snp_effect (@$effect) {

#      print Dumper( $trans_effect );

      $$snp_effect{ transcript_id } ||= "";
      
      my @effect_line;
      push @effect_line, "$$snp_effect{ external_name }/$$snp_effect{ stable_id }" if ($$snp_effect{ external_name } && 
										       $$snp_effect{ stable_id } );

      push @effect_line, "$$snp_effect{ stable_id }" if (!$$snp_effect{ external_name } && 
							 $$snp_effect{ stable_id } );

      push @effect_line, "" if (!$$snp_effect{ external_name } && 
							 !$$snp_effect{ stable_id } );


      push @effect_line,  "$$snp_effect{ transcript_id }";
      push @effect_line, $$snp_effect{ position } || "";
      push @effect_line, $$snp_effect{ cpos }     || "";
      push @effect_line, $$snp_effect{ ppos }     || "";
      push @effect_line, $$snp_effect{ rs_number } if ( $$snp_effect{ rs_number });
      
      print join("\t", @line, @effect_line) . "\n";
    }
  }
}



# 
# 
# 
# Kim Brugger (28 May 2010)
sub variation_effects {
  my ($var_features) = @_;


#  print "Nr of variants: " . @$var_features . "\n";

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

	  
	  if ( $$xref[0] ) {
	    $gene_res{ xref } = $$xref[0]->display_id;
	    $gene_res{ transcript_id} = join("/",$gene_res{ xref },$gene_res{ transcript_id});
	  }

	  $gene_res{ cpos } = "";
	  $gene_res{ ppos } = "";


	  $gene_res{ cpos } = "c.".$con->cdna_start if ( $con->cdna_start);

	  if ( $con->translation_start) {
	    my ( $old, $new ) = split("\/", $con->pep_allele_string);
	    
	    $new = $old if ( $string eq "SYNONYMOUS_CODING");
	    
	    $old = one2three( $old );
	    $new = one2three( $new );
	    

	    $gene_res{ ppos } = "p.$old".$con->translation_start . " $new";
	  }

	}

	$gene_res{ rs_number } = $existing_vf || "";
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
  if(defined $slice_hash{$chr}) {
    $slice = $slice_hash{$chr};
  }
 
  # if not create a new one
  else {
    
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
  my ( $chr, $SNP_pos, $ref, $alt) = @_;

  if ( ! $bam ) {
#    print STDERR "need a bam file for finding base distribution\n";
    return;
  }

  my %base_stats = ( A => 0, C => 0, G => 0, T => 0, N => 0);
  my %qual_stats;
  my $total = 0;
  
  open (my $st_pipe, "$samtools view $bam $chr:$SNP_pos-$SNP_pos | ") || die "Could not open samtools pipe: $!\n";

  while(<$st_pipe>) {
    chomp;
    my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");
    my @seq  = split("", $sequence);
    my @qual = split("", $quality);

    my $base_pos = $SNP_pos - $pos;
    my $base = $seq[ $base_pos ];
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


  my ($best, $score) = (1, 0);
  my %res;
  foreach my $base (sort {$base_stats{$b} <=> $base_stats{$a}} keys %base_stats ) {
    print "= sprintf(%.2f, $base_stats{$base}/$total*100); $chr $SNP_pos \n";
    my $perc = sprintf("%.2f", $base_stats{$base}/$total*100);
    my $qual = $qual_stats{$base} || 0;
    $res{$base} = "$base: $base_stats{$base}($perc%)/$qual";

    $score =  int($base_stats{$base} * $qual) if ( $base eq $alt);
    $score -= int($base_stats{$base} * $qual) if ( $base ne $alt && $base ne $ref);
  }

  $res{total} = $total;
  $res{score} = $score;

  
  return \%res;
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
# Kim Brugger (28 Apr 2010)
sub readin_cvf {
  my ($file, $min_score) = @_;
  open (my $in, $file) || die "Could not open '$file': $!\n";

  $min_score ||= 0;

  while(<$in>) {
    next if (/^\#/);
    
    my ($chr, $pos, $id, $ref_base, $alt_base, $mapq, $filter, $info) = split("\t");

#    print "$_";
    $alt_base = subtract_reference($alt_base, $ref_base);

    my %info_hash;
    foreach my $entry (split("\;", $info )) {
      my @f=split("\=", $entry); 
      $info_hash{$f[0]} = $f[1];
    }
    

    my $depth = $info_hash{ DP };

    next if ( $mapq < $min_score );
    
    $SNPs{$chr}{$pos}{GATK} = { depth        => $depth,
				mapping_qual => $mapq,
				alt_base     => $alt_base,
				pos => $pos,};

    $SNPs{$chr}{$pos}{ref_base} = $ref_base;
  }

}
 
# 
# 
# 
# Kim Brugger (28 Apr 2010)
sub readin_pileup {
  my ($file, $min_SNP_qual) = @_;
  open (my $in, $file) || die "Could not open '$file': $!\n";

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
 


