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

use lib '/usr/local/lib/ensembl-variation/modules/';
use lib '/usr/local/lib/ensembl-functgenomics/modules/';
use lib '/usr/local/lib/ensembl/modules/';
use lib '/usr/local/lib/bioperl/';

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use Scalar::Util qw/weaken/;

my %opts;
getopts('s:v:pb:Tq:m:HfrhnI', \%opts);
usage() if ( $opts{h});

my $species     = $opts{s} || "human";
my $buffer_size = 500;
my $host        = 'ensembldb.ensembl.org';
my $user        = 'anonymous';

# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user, -NO_CACHE => 0);
# get variation adaptors
my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
my $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');
my $sa  = $reg->get_adaptor($species, 'core', 'slice');
my $ga  = $reg->get_adaptor($species, 'core', 'gene');
my $afa = $reg->get_adaptor($species, 'funcgen', 'AnnotatedFeature');


my $weaken = 0;

my %slice_hash = ();


my $vcf         = $opts{v} || usage();
my $bam         = $opts{b};
my $regulation  = $opts{r} || 0;
my $pfam        = $opts{p} || 0;
my $from_36     = $opts{T} || 0;
my $min_mapq    = $opts{m} || undef;
my $no_ensembl  = $opts{n} || 0;
my $min_qual    = $opts{q} || undef;
my $full_report = $opts{f} || 0;
my $html_out    = $opts{H} || 0;
my $igv_links   = $opts{I} || 0;

my $dbsnp_link     = 'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=';
my $ens_snp_link   = 'http://www.ensembl.org/Homo_sapiens/Variation/Summary?v=';
my $ens_gene_link  = 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=';
my $ens_trans_link = 'http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=';

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

#readin_pileup( $pileup ) if ( $pileup);
readin_vcf( $vcf )      if ( $vcf);

my ($gatk_only, $pileup_only, $both_agree, $both_disagree) = (0, 0, 0, 0);
my $counter = 0;
foreach my $chr ( sort {$a cmp $b}  keys %SNPs ) {
  
  next if ($chr ne '1');


  my %res;
  my @vfs;

  foreach my $pos ( sort { $a <=> $b} keys %{$SNPs{$chr}} ) {
        
    
    my $position = "$chr:$pos";
    
    my @keys = keys %{$SNPs{$chr}{$pos}};

    @keys = grep(!/ref_base/, @keys);
    
    $res{$position}{callers} = join("/", @keys);

    if (1 || keys %{{ map {$SNPs{$chr}{$pos}{$_}{alt_base}, 1} @keys }} == 1) {
	
      my $key = "GATK";
      $res{$position}{ref_base}  = $SNPs{$chr}{$pos}{ref_base};
      $res{$position}{alt_base}  = $SNPs{$chr}{$pos}{$key}{alt_base};
      $res{$position}{qual}      = $SNPs{$chr}{$pos}{$key}{qual};
      $res{$position}{base_dist} = base_dist( $chr, $pos, $res{$position}{ref_base}, $res{$position}{alt_base});

      if ( ! $no_ensembl ) {

	my $Echr = $chr;
	$Echr =~ s/chr//;
	my $Epos = $pos;
	($Echr, $Epos) = remap($Echr, $pos, $pos) if ( $from_36 );
	
	next if ( ! $Echr );
	my $slice;
	$slice = $sa->fetch_by_region('chromosome', $Echr); 
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
	
	
	if ( 1 || @vfs >= $buffer_size) {
	  my $effects = variation_effects(\@vfs);
	  print_results( \%res, $effects);
	  undef @vfs;# = ();	  
	  undef %res;# = ();


	}
      }
    }
  }

  if ( @vfs || keys %res ) {
    my $effects = variation_effects(\@vfs) if ( ! $no_ensembl);
    print_results( \%res, $effects);
    @vfs = ();
    %res = ();
  }

}


print table_end() if ($html_out);



# 
# 
# 
# Kim Brugger (08 Jul 2010)
sub print_results {
  my ( $mapping, $effects ) = @_;
  
  if ( $full_report ) {
    print_fullreport($mapping, $effects );
  }    
  else { #if ( !$full_report ) {
    print_oneliner($mapping, $effects );
  }


}




#
# Creates a simple table, the function expects an array of arrays, and a border or not flag.
# 
# Kim Brugger (20 Oct 2003)
sub html_table {
  my ($cells,     # the cells as an array of arrays of values.
      $border,    # border or not.
      $padding,   # how big the cells are (padded around the text).
      $spacing,   # how the cells should be spaced
      $bgcolour,  # the colour of the table.
      $tablewidth, # how wide the table should be, this is a string, so we can handle both pixel width and percentages
      ) = @_;


  my $width = 0;
  foreach my $row (@$cells) {
    $width = @$row if ($width < @$row);
  }

#  my $return_string = table_start($border, $padding, $spacing, $bgcolour, $tablewidth);

  my $return_string = "";

  foreach my $row (@$cells) {
    
    $return_string .= "  <TR>";
    for (my $i=0; $i<$width;$i++) {
      $$row[$i] = "&nbsp;" if (!$$row[$i]);
      $return_string .= "<TD>$$row[$i]</TD>" 
    }
    $return_string .= "</TR>\n";
  }


  return $return_string;
}



#
# Creates an advanced table, the function expects an array of arrays, and a border or not flag.
# 
# Kim Brugger (20 Oct 2003)
sub advanced_table {
  my ($cells,      # the cells as an array of arrays of hashes, named according to the language:
                   #              %cell = (value, bgcolor, width, colspan, 
                   #                       rowspan, cellhalign, cellvalign);
      $border,     # border or not.
      $padding,    # how big the cells are (padded around the text).
      $spacing,    # how the cells should be spaced
      $bgcolour,   # the colour of the table.
      $spanning,   # Wether there are spanning cells, or if they should be "padded".
      $tablewidth, # how wide the table should be, this is a string, 
                   # so we can handle both pixel width and percentages
      $class       # Class so CSS can define the behaviour
      ) = @_;

  my $width = 0;
  foreach my $row (@$cells) {
    $width = @$row if ($width < @$row);
  }
  
#  my $return_string = table_start($border, $padding, $spacing, $bgcolour, $tablewidth, $class);
  my $return_string;
  foreach my $row (@$cells) {
    
    $return_string .= "  <TR>";
    for (my $i=0; $i<$width;$i++) {
      
      if ($$row[$i]) {
        $return_string .= "<TD";

        foreach my $key (keys %{$$row[$i]}) {
          next if ($key eq 'value');
          $return_string .= " $key='$$row[$i]{$key}'";
        }       
        $return_string .= ">$$row[$i]{'value'}</TD>" if ($$row[$i]{'value'});
      }
      else {
        $return_string .= "<TD>&nbsp;</TD>" if (!$spanning);
      }

    }
    $return_string .= "</TR>\n";
  }


  return $return_string;
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


sub table_start {
  my ($border, $padding, $spacing, $bgcolour, $width, $class) = @_;
  
  my $return_string .= "<TABLE";
  $return_string .= " border='$border'"        if $border;
  $return_string .= " cellspacing='$spacing'"  if $spacing;
  $return_string .= " cellpadding='$padding'"  if $padding;
  $return_string .= " bgcolor='$bgcolour'"     if $bgcolour;
  $return_string .= " width='$width'"          if $width;
  $return_string .= " class='$class'"          if $class;
  $return_string .= " >\n";

  return $return_string;
}

sub table_end {
  return "</TABLE>\n";

}



my ($printed_header) = (0);

# 
# 
# 
# Kim Brugger (08 Jul 2010)
sub print_oneliner {
  my ( $mapping, $effects ) = @_;

  my @res;
  if ( ! $printed_header++ ) {
    print table_start(1) if ( $html_out);
    
    my @annotations = ('position', 'change', 'score');
    push @annotations, ('depth', '', '', '', '', '') if ( $bam );
      
    push @annotations, ('gene', 'transcript', 'region', 'codon pos', 'AA change', 'external ref') if ( ! $no_ensembl);

    push @annotations, 'regulation' if ( $regulation );

    push @annotations, 'pfam' if ( $pfam );

    push @res, [@annotations];
  }

  if ( ! $effects ) {

    foreach my $name (sort keys %$mapping) {
      
      my @line;
      

      if ( $html_out && $igv_links) {

	if ($name =~ /^\d+\z/) {
	  push @line, "<a href='localhost:60151/locus=chr$name'> $name</a>";
	}
	else {
	  push @line,  "<a href='localhost:60151/locus=$name'> $name</a>";
	}
      }
      else {
	push @line, "$name";
      }

      push @line, "$$mapping{$name}{ref_base}>$$mapping{$name}{alt_base}";
      push @line, $$mapping{$name}{qual};
      if ( $$mapping{$name}{base_dist} ) {
	push @line, $$mapping{$name}{base_dist}{total};      
	map { push @line, $$mapping{$name}{base_dist}{$_} if ($$mapping{$name}{base_dist}{$_})} ( 'A', 'C', 'G', 'T', 'N');
      }
      
      push @res, \@line;
    }

  }
  else {

    foreach my $effect ( @$effects ) {
      
      my $name = $$effect[0]{ name };

      my @line;
      if ( $html_out && $igv_links) {

	if ($name =~ /^\d+:/) {
	  push @line, "<a href='http://localhost:60151/goto?locus=chr$name'> $name</a>";
	}
	else {
	  push @line,  "<a href='http://localhost:60151/goto?locus=$name'> $name</a>";
	}
      }
      else {
	push @line, "$name";
      }

      push @line, "$$mapping{$name}{ref_base}>$$mapping{$name}{alt_base}";
      push @line, $$mapping{$name}{qual};
      if ( $$mapping{$name}{base_dist} ) {
	push @line, $$mapping{$name}{base_dist}{total};      
	
	map { push @line, $$mapping{$name}{base_dist}{$_} if ($$mapping{$name}{base_dist}{$_})} ( 'A', 'C', 'G', 'T', 'N');
      }
      
      foreach my $snp_effect (@$effect) {
	
	my @effect_line;
	
	my $gene_id = "";
	
	$gene_id = "$$snp_effect{ external_name }/$$snp_effect{ stable_id }" if ($$snp_effect{ external_name } && 
										 $$snp_effect{ stable_id } );
	
	$gene_id = "$$snp_effect{ stable_id }" if (!$$snp_effect{ external_name } && 
						   $$snp_effect{ stable_id } );
	
	$gene_id = "<a href='$ens_gene_link$$snp_effect{ stable_id }'>$gene_id</a>"
	    if ( $gene_id && $html_out && $$snp_effect{ stable_id });
	
	push @effect_line, $gene_id;
	
	
	my $trans_id = "";
	
	$trans_id = "$$snp_effect{ xref }/$$snp_effect{ transcript_id }" if ($$snp_effect{ xref } && 
									     $$snp_effect{ transcript_id } );
	
	$trans_id = "$$snp_effect{ transcript_id }" if (!$$snp_effect{ xref } && 
							$$snp_effect{ transcript_id } );
	
	$trans_id = "<a href='$ens_trans_link$$snp_effect{ transcript_id }'>$trans_id</a>"
	    if ( $trans_id && $html_out && $$snp_effect{ transcript_id });
	
	
	push @effect_line, $trans_id;
	push @effect_line, $$snp_effect{ position } || "";
	push @effect_line, $$snp_effect{ cpos }     || "";
	push @effect_line, $$snp_effect{ ppos }     || "";
	

	$$snp_effect{ rs_number } = "<a href='$dbsnp_link$$snp_effect{ rs_number }'>$$snp_effect{ rs_number }</a>"
	    if ( $$snp_effect{ rs_number } && $html_out && $$snp_effect{ rs_number } =~ /rs\d+/ && $$snp_effect{ rs_number } !~ /href/);

	$$snp_effect{ rs_number } = "<a href='$ens_snp_link$$snp_effect{ rs_number }'>$$snp_effect{ rs_number }</a>"
	    if ( $$snp_effect{ rs_number } && $html_out && $$snp_effect{ rs_number } =~ /ENSSNP\d+/ && $$snp_effect{ rs_number } !~ /href/);

      push @effect_line, $$snp_effect{ rs_number } || "";

      push @effect_line, $$snp_effect{ pfam } || "" if ( $pfam);
      push @effect_line, $$snp_effect{ regulation } || "" if ( $regulation);

	push @res, [@line, @effect_line];
      }
    }
  }

  if (  $html_out ) {
    print html_table(\@res, 1);
  }
  else {
    print text_table(\@res, 1);
  }
  
}



# 
# 
# 
# Kim Brugger (08 Jul 2010)
sub print_fullreport {
  my ( $mapping, $effects ) = @_;

  my @res;
  if ( ! $printed_header++ ) {
    print table_start(1) if ( $html_out);
    
#    my @annotations = ('position', 'change', 'score');
#    push @annotations, ('depth', '', '', '', '', '') if ( $bam );
      
#    push @annotations, ('gene', 'transcript', 'region', 'codon pos', 'AA change', 'external ref');

#    push @annotations, ('pfam', 'regulation') if ( $full );

#    push @res, [@annotations];
  }

  foreach my $effect ( @$effects ) {

    my $name = $$effect[0]{ name };

    my @line;

    if ( $html_out && $igv_links) {

      if ($name =~ /^\d+:/) {
	push @line, "<a href='http://localhost:60151/goto?locus=chr$name'> $name</a>";
      }
      else {
	push @line,  "<a href='http://localhost:60151/goto?locus=$name'> $name</a>";
      }
    }
    else {
      push @line, "$name";
    }

    push @line, "$$mapping{$name}{ref_base}>$$mapping{$name}{alt_base}";


    push @line, $$mapping{$name}{qual};
    if ( $$mapping{$name}{base_dist} ) {
      push @line, $$mapping{$name}{base_dist}{total};      
    
      map { push @line, $$mapping{$name}{base_dist}{$_} if ($$mapping{$name}{base_dist}{$_})} ( 'A', 'C', 'G', 'T', 'N');
    }

    push @line, $$effect[0]{ rs_number } || "";
    push @line, $$effect[0]{ regulation } || "" if ( $regulation );
    push @res, \@line;
    
    foreach my $snp_effect (@$effect) {

      my @effect_line = ("","");


      my $gene_id = "";

      $gene_id = "$$snp_effect{ external_name }/$$snp_effect{ stable_id }" if ($$snp_effect{ external_name } && 
									       $$snp_effect{ stable_id } );

      $gene_id = "$$snp_effect{ stable_id }" if (!$$snp_effect{ external_name } && 
						 $$snp_effect{ stable_id } );

      $gene_id = "<a href='$ens_gene_link$$snp_effect{ stable_id }'>$gene_id</a>"
	  if ( $gene_id && $html_out && $$snp_effect{ stable_id });

      push @effect_line, $gene_id;


      my $trans_id = "";

      $trans_id = "$$snp_effect{ xref }/$$snp_effect{ transcript_id }" if ($$snp_effect{ xref } && 
									   $$snp_effect{ transcript_id } );

      $trans_id = "$$snp_effect{ transcript_id }" if (!$$snp_effect{ xref } && 
						      $$snp_effect{ transcript_id } );

      $trans_id = "<a href='$ens_trans_link$$snp_effect{ transcript_id }'>$trans_id</a>"
	  if ( $trans_id && $html_out && $$snp_effect{ transcript_id });


      push @effect_line, $trans_id;
      push @effect_line, $$snp_effect{ position } || "";
      push @effect_line, $$snp_effect{ cpos }     || "";
      push @effect_line, $$snp_effect{ ppos }     || "";
      

      $$snp_effect{ rs_number } = "<a href='$dbsnp_link$$snp_effect{ rs_number }'>$$snp_effect{ rs_number }</a>"
	  if ( $$snp_effect{ rs_number } && $html_out && $$snp_effect{ rs_number } =~ /rs\d+/ && $$snp_effect{ rs_number } !~ /href/);

      $$snp_effect{ rs_number } = "<a href='$ens_snp_link$$snp_effect{ rs_number }'>$$snp_effect{ rs_number }</a>"
	  if ( $$snp_effect{ rs_number } && $html_out && $$snp_effect{ rs_number } =~ /ENSSNP\d+/ && $$snp_effect{ rs_number } !~ /href/);


      push @effect_line, $$snp_effect{ pfam } || "";
      push @res, [@effect_line];
    }
  }
    

  if (  $html_out ) {
    print html_table(\@res, 1);
  }
  else {
    print text_table(\@res, 1);
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
    my $regulations;
    
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

      if ( $regulation ) {
	
	my $features = $afa->fetch_all_by_Slice($fs); 
	if ( @$features ) {
	  my @types;
	  foreach my $f ( @{$features}) {
	    push @types, $f->display_label();
	  }
	  
	  my %saw;
	  $regulations = join("/", grep(!$saw{$_}++,@types));
	}
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

	  
	  $gene_res{ xref } = $$xref[0]->display_id
	      if ( $$xref[0] );

	  $gene_res{ cpos } = "";
	  $gene_res{ ppos } = "";


	  $gene_res{ cpos } = "c.".$con->cdna_start if ( $con->cdna_start);

	  if ( $con->translation_start) {
	    my ( $old, $new ) = ("","");
	    ( $old, $new ) = split("\/", $con->pep_allele_string);
	    
	    $new = $old if ( !$new || $new eq "");
	    $old = one2three( $old );
	    $new = one2three( $new );
	    $gene_res{ ppos } = "p.$old".$con->translation_start . " $new";

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
	    weaken($protein) if ($weaken);
	    weaken($prot_feats) if ($weaken);
	  }
	  weaken($gene) if ($weaken);
	  weaken($xref) if ($weaken);
	}

	$gene_res{ regulation } = $regulations; 
	$gene_res{ rs_number }  = $existing_vf || "";
	push @{$res[$feature]}, \%gene_res;

      }
      weaken($con) if ($weaken);
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
  my ( $chr, $SNP_pos) = @_;

  if ( ! $bam ) {
#    print STDERR "need a bam file for finding base distribution\n";
    return;
  }

  my %base_stats = ( A => 0, C => 0, G => 0, T => 0, N => 0);
  my %qual_stats;
  my $total = 0;

  open (my $st_pipe, "$samtools view $bam $chr:$SNP_pos-$SNP_pos | ") || die "Could not open samtools pipe: $!";

  while(<$st_pipe>) {
    chomp;
    my ($read, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $sequence, $quality, $opts) = split("\t");
    next if ( defined $min_mapq && $min_mapq > $mapq);

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
  foreach my $base (sort {$base_stats{$b} <=> $base_stats{$a}} keys %base_stats ) {
    my $perc = sprintf("%.2f", $base_stats{$base}/$total*100);
    my $qual = $qual_stats{$base} || 0;
    $res{$base} = "$base: $base_stats{$base}($perc%)/$qual";

  }

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
sub readin_vcf {
  my ($file) = @_;
  open (my $in, $file) || die "Could not open '$file': $!";


  while(<$in>) {
    next if (/^\#/);
    
    my ($chr, $pos, $id, $ref_base, $alt_base, $qual, $filter, $info) = split("\t");

#    next if ( defined $min_qual && $min_qual > $qual);

    next if ( $filter ne "PASS" && $filter ne "." );

#    print "$_";
    $alt_base = subtract_reference($alt_base, $ref_base);

    my %info_hash;
    foreach my $entry (split("\;", $info )) {
      my @f=split("\=", $entry); 
      $info_hash{$f[0]} = $f[1];
    }
    

    my $depth = $info_hash{ DP };

    
    $SNPs{$chr}{$pos}{GATK} = { depth     => $depth,
				qual      => $qual,
				alt_base  => $alt_base,
				pos       => $pos,};

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
 



# 
# 
# 
# Kim Brugger (09 Jul 2010)
sub usage {

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

=item B<-r>:

Extracts regulation information from ensembl

=item B<-v F<vcf file>>: 

The infile that is to be converted into a nice report.

=item B<-T>: 

The mapping was done on NCBI36/hg18 and the coordinates should be transformed to GRCh37/hg19.

=back

=head1 NOTES

=over

=item B<-s>: will change the species, but changing that will probably break about just everything.


=cut
