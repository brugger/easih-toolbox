#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (28 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


use lib '/software/lib/ensembl-variation/modules/';
use lib '/software/lib/ensembl-functgenomics/modules/';
use lib '/software/lib/ensembl/modules/';
use lib '/software/lib/bioperl/';

use strict;
use Getopt::Std;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;

my %opts;
getopts('b:d:TfHh', \%opts);
usage() if ( $opts{h});

my $species     = "human";
my $buffer_size = 500;
my $host        = 'ensembldb.ensembl.org';
my $user        = 'anonymous';

my $host        = 'localhost';
$host           = 'mgpc17';
my $user        = 'easih_ro';


# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);
# get variation adaptors
my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
my $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');
my $sa = $reg->get_adaptor($species, 'core', 'slice');
my $ga = $reg->get_adaptor($species, 'core', 'gene');

my $input       = $opts{b} || usage();
my $from_36     = $opts{T} || 0;
my $full_report = $opts{f} || 0;
my $html_out    = $opts{H} || 0;

my $ens_gene_link = 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=';
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


my $dips_n_gaps = readin_input( $input );
foreach my $chr ( sort {$a cmp $b}  keys %$dips_n_gaps ) {
  

  foreach my $start ( sort { $a <=> $b} keys %{$$dips_n_gaps{$chr}} ) {

    my $dip = $$dips_n_gaps{$chr}{ $start };
    my $end   = $$dip{end};
    my $position = "$chr:$start-$end";

    my @line;
    push @line, $position, $$dip{depth};
    my $effects = region_effects($chr, $start, $end);


    print_results( \@line, $effects );
    
  }


#  last;
}


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

my ($printed_header) = (0);


# 
# 
# 
# Kim Brugger (05 Aug 2010)
sub print_fullreport {
  my ( $mapping, $effects ) = @_;

  if ( ! $printed_header++ ) {
    print table_start(1) if ( $html_out);
  }
  
  my @res;

  push @res, [@$mapping]; 

  if ( @$effects ) {
      
    foreach my $effect ( @$effects ) {
      
      my @effect_line;
      my $gene_id;
      $gene_id = "$$effect{ external_name }/$$effect{ stable_id }" if ($$effect{ external_name } && 
								       $$effect{ stable_id } );

      $gene_id = "$$effect{ stable_id }" if (!$$effect{ external_name } && 
					     $$effect{ stable_id } );
	


      $gene_id = "$$effect{ external_name }" if ($$effect{ external_name });

      $gene_id = "<a href='$ens_gene_link$$effect{ stable_id }'>$gene_id</a>"
	  if ( $gene_id && $html_out && $$effect{ stable_id });
	
      push @effect_line, $gene_id;

      my $trans_id = "$$effect{ xref }/$$effect{ transcript_id }" if ($$effect{ xref } && 
								   $$effect{ transcript_id } );
      
      $trans_id = "$$effect{ transcript_id }" if (!$$effect{ xref } && 
						  $$effect{ transcript_id } );
      
      $trans_id = "<a href='$ens_trans_link$$effect{ transcript_id }'>$trans_id</a>"
	  if ( $trans_id && $html_out && $$effect{ transcript_id });
	
      push @effect_line, $trans_id;
#      push @effect_line, ($$effect{ position } || "");
#      push @effect_line, ($$effect{ cpos } || "");
#      push @effect_line, ($$effect{ ppos } || "");
      
      push @res, ["","", @effect_line];
    }
  }

  if (  $html_out ) {
    print html_table(\@res, 1);
  }
  else {
    print text_table(\@res, 1);
  }

}




sub print_oneliner {
  my ( $mapping, $effects ) = @_;

  my @res;
  if ( ! $printed_header++ ) {
    print table_start(1) if ( $html_out);
    
#    push @res, ['position', 'depth', 'gene', 'transcript', 'region', 'codon pos', 'protein pos' ];
    push @res, ['position', 'depth', 'gene' ];

  }

  if ( @$effects ) {
      
    foreach my $effect ( @$effects ) {
      
      my @effect_line;
#      push @effect_line, "$$effect{ external_name }/$$effect{ stable_id }", "$$effect{ transcript_id }";
      push @effect_line, "$$effect{ external_name }";
      push @effect_line, ($$effect{ position } || "");
      push @effect_line, ($$effect{ cpos } || "");
      push @effect_line, ($$effect{ ppos } || "");
      
      push @res, [ @$mapping, @effect_line];
    }
  }
  else {
    push @res, [ @$mapping];
  }


  if (  $html_out ) {
    print html_table(\@res, 1);
  }
  else {
    print text_table(\@res, 1);
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
# 
# 
# Kim Brugger (28 May 2010)
sub region_effects {
  my ( $chr, $start, $end ) = @_;

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
    -allele_string => "",
    -strand => $strand,
    -map_weight => 1,
    -adaptor => $vfa,           # we must attach a variation feature adaptor
    -variation_name => $chr.'_'.$start.'_'.$end,
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

	  
	  $gene_res{ xref } = $$xref[0]->display_id
	      if ( $$xref[0] );

	  $gene_res{ cpos } = "";
	  $gene_res{ ppos } = "";


	  $gene_res{ position } = $string;
	  $gene_res{ cpos } = "c.".$con->cdna_start if ( $con->cdna_start);

	  if ( $con->translation_start && $con->pep_allele_string) {
	    my ( $old, $new ) = split("\/", $con->pep_allele_string);
	    
#	    print " --> " . $con->pep_allele_string . "\n";

	    $new ||= $old;

	    $old = one2three( $old );
	    $new = one2three( $new );
	    

	    $gene_res{ ppos } = "p.$old".$con->translation_start;
	    $gene_res{ position } = "WITHIN_CODING_REGION";
	  }

	  push @res, \%gene_res;
	} 
      }
      last;
    }
  }

#  print Dumper( \@res );
  
  return \@res;
}





# 
# 
# 
# Kim Brugger (02 Jun 2010)
sub one2three {
  my ( $aminoacids) = @_;
  
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

  return join(",", map{ $trans{ $_ }} split("", $aminoacids));

#  return $trans{ $aminoacid } if ($trans{ $aminoacid });
#  return $aminoacid;
}





# 
# 
# 
# Kim Brugger (28 Apr 2010)
sub readin_input {
  my ($file) = @_;

  my %input = ();
  open (my $in, $file) || die "Could not open '$file': $!\n";
  while(<$in>) {

    next if (/^\#/);
    chomp;
    
    my ($region, $depth) = split("\t");

    my ($name, $start, $end) = $region =~ /^(\w+?):(\d+)-(\d+)/;
    $input{ $name }{ $start } = { end    => $end,
				  depth  => $depth};
    
  }



  return \%input;
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

indel_report turns a vcf file into a nice report, annotated with Ensembl gene information 

=head1 OPTIONS


=over

=item B<-b F<bed file>>: 

The bedfile resulting from the gap/dip calling.

=item B<-f>: 

Prints out a multi-line report, otherwise it is done on a oneline basis (good for excel)

=item B<-H>: 

Prints a HTML report, default is a tab-separated one.

=item B<-T>: 

The mapping was done on NCBI36/hg18 and the coordinates should be transformed to GRCh37/hg19.

=back



=cut
