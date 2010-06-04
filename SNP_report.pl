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

my $gatk        = shift;
my $pileup      = shift;
my $bam         = shift || "";
my $from_36     = 1;
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

my ($gatk_only, $pileup_only, $both_agree, $both_disagree) = (0,0,0, 0);
foreach my $chr ( sort {$a cmp $b}  keys %SNPs ) {
  
  $chr = 'chrX';

  foreach my $pos ( sort { $a <=> $b} keys %{$SNPs{$chr}} ) {

    my @fields;
    
    my $position = "$chr:$pos";
    
    push @fields, "chr:pos", "Ref base", "All. base", "SNP caller(s)", "Base stats #base(percent)/avg score", "\n$chr:$pos";

    push @fields, $SNPs{$chr}{$pos}{ref_base};

    my @keys = keys %{$SNPs{$chr}{$pos}};

    @keys = grep(!/ref_base/, @keys);
    
    $SNPs{$chr}{$pos}{callers} = join("/", @keys);

    if (keys %{{ map {$SNPs{$chr}{$pos}{$_}{alt_base}, 1} @keys }} == 1) {
	
      my $key = $keys[0];

      $SNPs{$chr}{$pos}{alt_base} = $SNPs{$chr}{$pos}{$key}{alt_base};

      $SNPs{$chr}{$pos}{base_dist}  = base_dist( $chr, $pos, $SNPs{$chr}{$pos}{ref_base}, $SNPs{$chr}{$pos}{alt_base});
      $SNPs{$chr}{$pos}{snp_effect} = snp_effect($chr, $pos, $pos, "$SNPs{$chr}{$pos}{ref_base}/$SNPs{$chr}{$pos}{alt_base}");
    }
    else {

      my @snps;     
      foreach my $key ( @keys ) {
	push @snps, $SNPs{$chr}{$pos}{$key}{alt_base};
      }
      
      $SNPs{$chr}{$pos}{alt_base}   = join("/", @snps);
      $SNPs{$chr}{$pos}{base_dist}  = base_dist( $chr, $pos);;

      foreach my $snp ( @snps ) {
#	push @{$SNPs{$chr}{$pos}{snp_effects}}, snp_effect($chr, $pos, $pos, $snp );
      }
    }

#    print  Dumper( $SNPs{$chr}{$pos} );

    if ( ! $full_report && ! $html_out ) {
      
      my @line;
      push @line, "$chr:$pos", "$SNPs{$chr}{$pos}{ref_base}>$SNPs{$chr}{$pos}{alt_base}";
      push @line, ${$SNPs{$chr}{$pos}}{base_dist}{score};
      push @line, ${$SNPs{$chr}{$pos}}{base_dist}{total};      

      map { push @line, ${$SNPs{$chr}{$pos}}{base_dist}{$_} if (${$SNPs{$chr}{$pos}}{base_dist}{$_})} ( 'A', 'C', 'G', 'T', 'N');
      push @line, ${$SNPs{$chr}{$pos}}{callers};

      if ($SNPs{$chr}{$pos}{snp_effect} ) {
      
	foreach my $snp_effect ( @{$SNPs{$chr}{$pos}{snp_effect}} ) {
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
#    exit;
    

  }

  last;
}




# 
# 
# 
# Kim Brugger (29 Apr 2010)
sub base_dist {
  my ( $chr, $SNP_pos, $ref, $alt) = @_;

  if ( ! $bam ) {
    print STDERR "need a bam file for finding base distribution\n";
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

  my ($first, $second, $first_qual, $second_qual) = (0,0, 0, 0);
  foreach my $key ( keys %qual_stats ) {
    
    my $sum = eval join '+', @{$qual_stats{$key}};
    my $count = @{$qual_stats{$key}};
    my $avg_qual      =  sprintf("%.2f",$sum / $count);
    $avg_qual ||= 0;


    if ( $count ) {
      if ( ! $first || $first < $count ) { 
      
	if ( $first ) {
	  $second      = $first;
	  $second_qual = $first_qual;
	}
	$first      = $count;
	$first_qual = $avg_qual;
      }
      elsif (! $second || $second < $count ) {
	$second      = $count;
	$second_qual = $avg_qual;
      }      
    }
  
    $qual_stats{$key} = $avg_qual;
  }



  my $ratio = $second/( $first + $second) * 100;
  my $ignore_second = 0;
  if ( $ratio < 10 ) {    
    $ignore_second++;
  }
  elsif ( $ratio > 15 && $ratio < 35 ) {
    $ignore_second++;    
  }
  elsif ( $ratio > 40 && $ratio < 60 ) {
    $ignore_second++;    
  }
  elsif ($ratio > 60 && $ratio < 83 ) {
    $ignore_second++;    
  }

  my ($best, $score) = (1, 0);
  my %res;
  foreach my $base (sort {$base_stats{$b} <=> $base_stats{$a}} keys %base_stats ) {
    my $perc = sprintf("%.2f", $base_stats{$base}/$total*100);
    my $qual = $qual_stats{$base} || 0;
    $res{$base} = "$base: $base_stats{$base}($perc%)/$qual";
    if ( $best ) {
      $score = int($base_stats{$base} * $qual);
      $best = 0;
    }
    elsif ( $ignore_second ) {
      $ignore_second = 0;
    }
    else {
      $score -= int($base_stats{$base} * $qual);
    }
  }

  $res{total} = $total;
  $res{score} = $score;

  
  return \%res;
}



# 
# 
# 
# Kim Brugger (28 May 2010)
sub snp_effect {
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
  
  return \@res;
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
sub readin_geli {
  my ($file, $min_depth) = @_;
  open (my $in, $file) || die "Could not open '$file': $!\n";

  $min_depth ||= 0;

  while(<$in>) {
    next if (/^\#/);
    
    my ($chr, $pos, $ref_base, $depth, $mapq, $genotype, $best_lod) = split(" ");

    next if ( $depth < $min_depth );
    
    $SNPs{$chr}{$pos}{GATK} = { mapping_qual => $mapq,
				genotype     => $genotype,
				lodscore     => $best_lod,
				pos          => $pos};
    
    $SNPs{$chr}{$pos}{ref_base} = $ref_base;

  }

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
 


