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

my $gatk   = shift;
my $pileup = shift;
my $bam    = shift || "";

my $from_36 = 1;
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
    push @fields, "chr:pos", "Ref base", "All. base", "SNP caller(s)", "Base stats #base(percent)/avg score", "\n$chr:$pos";

    push @fields, $SNPs{$chr}{$pos}{ref_base};

    my @keys = keys %{$SNPs{$chr}{$pos}};

    @keys = grep(!/ref_base/, @keys);

    if ( @keys == 1 ) {
      my $key = $keys[0];

      push @fields, $SNPs{$chr}{$pos}{$key}{alt_base};
      push @fields, $key;

      push @fields, base_dist( $chr, $pos);;

      push @fields, snp_effect($chr, $pos, $pos, $SNPs{$chr}{$pos}{$key}{alt_base});

    }
    else {

      my @snps;

      if (keys %{{ map {$SNPs{$chr}{$pos}{$_}{alt_base}, 1} @keys }} == 1) {
	
	my $key = $keys[0];
	
	push @snps, $SNPs{$chr}{$pos}{$key}{alt_base};

	push @fields, $snps[0];
	push @fields, join("/", @keys);

	push @fields, base_dist( $chr, $pos);;
	push @fields, snp_effect($chr, $pos, $pos, $snps[0] );


      }
      else {
	foreach my $key ( @keys ) {
	  
	  push @snps, $SNPs{$chr}{$pos}{$key}{alt_base};
	}
      
	push @fields, join("/", @snps);
	push @fields, join("/", @keys);

	push @fields, base_dist( $chr, $pos);;
	push @fields, snp_effect($chr, $pos, $pos, $snps[0] );
      }


    }
    
    

    print join("\t", @fields) . "\n";

  }

  last;
}




# 
# 
# 
# Kim Brugger (29 Apr 2010)
sub base_dist {
  my ( $chr, $SNP_pos) = @_;

  if ( ! $bam ) {
    print STDERR "need a bam file for finding base distribution\n";
    return;
  }

  my %base_stats = ( A => 0, C => 0, G => 0, T => 0);
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
  elsif ( $ratio > 40 && $ratio < 60 ) {
    $ignore_second++;    
  }
  elsif ($ratio > 60 && $ratio < 80 ) {
    $ignore_second++;    
  }

  my ($best, $score) = (1, 0);
  my @res;
  foreach my $base (sort {$base_stats{$b} <=> $base_stats{$a} } keys %base_stats ) {
    my $perc = sprintf("%.2f", $base_stats{$base}/$total*100);
    my $qual = $qual_stats{$base} || 0;
    push @res, "$base: $base_stats{$base}($perc%)/$qual";
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


  push @res, "SNP score: $score";

  
  return @res;
  
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
  push @res, "\n\t";
  push @res, "Gene","Transcript", "Consequence","Position in cDNA","Position in protein", "Amino acid change","Corresponding Variation", "\n\t";
  
  my $added_entries = 0;

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
    my $existing_vf = "N/A";
    
    if(defined($new_vf->adaptor->db)) {
      my $fs = $new_vf->feature_Slice;
      if($fs->start > $fs->end) {
	($fs->{'start'}, $fs->{'end'}) = ($fs->{'end'}, $fs->{'start'});
      }
      foreach my $existing_vf_obj(@{$new_vf->adaptor->fetch_all_by_Slice($fs)}) {
	$existing_vf = $existing_vf_obj->variation_name
	    if $existing_vf_obj->seq_region_start == $new_vf->seq_region_start
	    and $existing_vf_obj->seq_region_end == $new_vf->seq_region_end;
      }
    }
		
    # the get_all_TranscriptVariations here now just retrieves the
    # objects that were attached above - it doesn't go off and do
    # the calculation again		
    foreach my $con (@{$new_vf->get_all_TranscriptVariations}) {
      foreach my $string (@{$con->consequence_type}) {
	
	next if ( $string eq "INTERGENIC");

	if($con->cdna_start && $con->cdna_end && $con->cdna_start > $con->cdna_end) {
	  ($con->{'cdna_start'}, $con->{'cdna_end'}) = ($con->{'cdna_end'}, $con->{'cdna_start'});
	}
	
	if($con->translation_start &&  $con->translation_end && $con->translation_start > $con->translation_end) {
	  ($con->{'translation_start'}, $con->{'translation_end'}) = ($con->{'translation_end'}, $con->{'translation_start'});
	}

	push @res, ($con->transcript ? $ga->fetch_by_transcript_stable_id($con->transcript->stable_id)->stable_id : "N/A");
	push @res, ($con->transcript ? $con->transcript->stable_id : "N/A");
	push @res, $string;
	push @res, ($con->cdna_start ? $con->cdna_start.($con->cdna_end eq $con->cdna_start ? "" : "-".$con->cdna_end) : "N/A");
	push @res, ($con->translation_start ? $con->translation_start.($con->translation_end eq $con->translation_start ? "" : "-".$con->translation_end) : "N/A");
	push @res, ($con->pep_allele_string ? $con->pep_allele_string : "N/A");
	push @res, $existing_vf;
	$added_entries++;
      }
      push @res, "\n\t";
    }
  }

#  print "\n" . join("\t", @res) . "\n";

  return @res if ( $added_entries );
  return "\n";

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

    next if ( $cons_SNP_qual < $min_SNP_qual);
    
    $SNPs{$chr}{$pos}{samtools} = { depth        => $depth,
				mapping_qual => $cons_SNP_qual,
				alt_base     => $alt_base,
				pos => $pos,};

    $SNPs{$chr}{$pos}{ref_base} = $ref_base;
  }
  
}
 


