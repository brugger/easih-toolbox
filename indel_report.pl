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

# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);
# get variation adaptors
my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
my $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');
my $sa = $reg->get_adaptor($species, 'core', 'slice');
my $ga = $reg->get_adaptor($species, 'core', 'gene');

my $bed         = $opts{b} || usage();
my $min_depth   = $opts{d} || 0;
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


my $indels = readin_bed( $bed )      if ( $bed);
foreach my $chr ( sort {$a cmp $b}  keys %$indels ) {
  

  foreach my $start ( sort { $a <=> $b} keys %{$$indels{$chr}} ) {

    my $indel = $$indels{$chr}{ $start };
    my $end   = $$indel{end};
    my $position = "$chr:$start-$end";

#    print Dumper( $indel );

    next if ( $$indel{depth} < 20);

    my @line;
    push @line, $position;
    push @line, $$indel{type};
    push @line, $$indel{variation};
    push @line, $$indel{support} . "/". $$indel{depth};
    
    my $effects = indel_effect($chr, $start, $end, "$$indel{variation}/$$indel{type}");


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
	
      $gene_id = "<a href='$ens_gene_link$$effect{ stable_id }'>$gene_id</a>"
	  if ( $gene_id && $html_out && $$effect{ stable_id });
	
      push @effect_line, $gene_id;

      my $trans_id = "";
	
      $trans_id = "$$effect{ xref }/$$effect{ transcript_id }" if ($$effect{ xref } && 
								   $$effect{ transcript_id } );
      
      $trans_id = "$$effect{ transcript_id }" if (!$$effect{ xref } && 
						  $$effect{ transcript_id } );
      
      $trans_id = "<a href='$ens_trans_link$$effect{ transcript_id }'>$trans_id</a>"
	  if ( $trans_id && $html_out && $$effect{ transcript_id });
	
	
      push @effect_line, $trans_id;
      push @effect_line, ($$effect{ position } || "");
      push @effect_line, ($$effect{ cpos } || "");
      push @effect_line, ($$effect{ ppos } || "");
      
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
    
    push @res, ['position', 'change', 'base(s)', 'evidence/depth','gene', 'transcript', 'region', 'codon pos', 'protein pos' ];

  }

  if ( @$effects ) {
      
    foreach my $effect ( @$effects ) {
      
      my @effect_line;
      push @effect_line, "$$effect{ external_name }/$$effect{ stable_id }", "$$effect{ transcript_id }";
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
#  $allele_string = "C/GGGG";
#  print "AS :: $allele_string\n";

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

	  
	  $gene_res{ xref } = $$xref[0]->display_id
	      if ( $$xref[0] );

	  $gene_res{ cpos } = "";
	  $gene_res{ ppos } = "";


	  $gene_res{ position } = $string;
	  $gene_res{ cpos } = "c.".$con->cdna_start if ( $con->cdna_start);

	  if ( $con->translation_start && $con->pep_allele_string) {
	    my ( $old, $new ) = split("\/", $con->pep_allele_string);
	    
	    print " --> " . $con->pep_allele_string . "\n";

	    $new ||= $old;

	    $old = one2three( $old );
	    $new = one2three( $new );
	    

	    $gene_res{ ppos } = "p.$old".$con->translation_start . " $new";
	  }

	  push @res, \%gene_res;
	} 
      }
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
				  depth     => $depth,
				  support   => $support,
				  change    => $change };
    
  }



  return \%indels;
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

The bedfile that the SNP calling was based on.

=item B<-f>: 

Prints out a multi-line report, otherwise it is done on a oneline basis (good for excel)

=item B<-H>: 

Prints a HTML report, default is a tab-separated one.

=item B<-d>: 

Minumum depth required when reporting an indel

=item B<-T>: 

The mapping was done on NCBI36/hg18 and the coordinates should be transformed to GRCh37/hg19.

=back



=cut
