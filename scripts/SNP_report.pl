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


use lib '/home/easih/lib/ensembl-variation/modules/';
use lib '/home/easih/lib/ensembl-functgenomics/modules/';
use lib '/home/easih/lib/ensembl/modules/';
use lib '/home/easih/lib/bioperl/';


use lib '/home/kb468/easih-toolbox/modules/';
use EASIH::SNPs;


use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;


my %opts;
getopts('s:v:pb:Tq:m:HfrhnILEg', \%opts);
perldoc() if ( $opts{h});

my $species     = $opts{s} || "human";
my $buffer_size = 100;
my $host        = 'ensembldb.ensembl.org';
my $user        = 'anonymous';

my $exit_count = 10;

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
my $at_ENSEMBL  = $opts{E} || 0;
my $grantham    = $opts{g} || 0;

# get registry
my $reg = 'Bio::EnsEMBL::Registry';
if ( $at_ENSEMBL ) {
  $reg->load_registry_from_db(-host => $host,-user => $user, -NO_CACHE => 0,);
}
else {
  $reg->load_registry_from_db(-host => "localhost",-user => "easih_ro", -NO_CACHE => 0);
}
# get variation adaptors
my $vfa = $reg->get_adaptor($species, 'variation', 'variationfeature');
my $tva = $reg->get_adaptor($species, 'variation', 'transcriptvariation');
my $sa  = $reg->get_adaptor($species, 'core', 'slice');
my $ga  = $reg->get_adaptor($species, 'core', 'gene');
my $afa = $reg->get_adaptor($species, 'funcgen', 'AnnotatedFeature');




my %slice_hash = ();
my ($sth_dbsnp, $sth_pop);
my $use_local_dbsnp = $opts{L} || 0;
if ($use_local_dbsnp) {
  require DBI;
  $use_local_dbsnp = DBI->connect('DBI:mysql:dbsnp132', 'easih_ro') || die "Could not connect to database: $DBI::errstr";
  $sth_dbsnp = $use_local_dbsnp->prepare("SELECT rs, flags FROM snp WHERE chr=? AND pos = ?");
  $sth_pop   = $use_local_dbsnp->prepare("select * from population where rs=?;");

}


my $hgmd     = DBI->connect('DBI:mysql:hgmd', 'easih_ro') || print STDERR "Could not connect to database hgmd: $DBI::errstr";
my $sth_hgmd = $hgmd->prepare("SELECT rs, chr, pos FROM snp WHERE (chr=? AND pos = ?) OR rs=?;");


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

my %grch37_remapping = ();

my ($gatk_only, $pileup_only, $both_agree, $both_disagree) = (0, 0, 0, 0);
my $counter = 10;
foreach my $chr ( sort {$a cmp $b}  keys %SNPs ) {
  

  my %res;
  my @vfs;

  foreach my $pos ( sort { $a <=> $b} keys %{$SNPs{$chr}} ) {
        
    
    my $position = "$chr:$pos";
    
    my @keys = keys %{$SNPs{$chr}{$pos}};

    @keys = grep(!/ref_base/, @keys);
    
    $res{$position}{callers} = join("/", @keys);

    if ( keys %{{ map {$SNPs{$chr}{$pos}{$_}{alt_base}, 1} @keys }} == 1) {
	
      my $key = "GATK";
      $res{$position}{ref_base}  = $SNPs{$chr}{$pos}{ref_base};
      $res{$position}{alt_base}  = $SNPs{$chr}{$pos}{$key}{alt_base};
      $res{$position}{qual}      = $SNPs{$chr}{$pos}{$key}{qual};
      $res{$position}{filter}    = $SNPs{$chr}{$pos}{$key}{filter};
      $res{$position}{base_dist} = base_dist( $chr, $pos, $res{$position}{ref_base}, $res{$position}{alt_base});

      if ( ! $no_ensembl ) {

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
	  my $effects = variation_effects(\@vfs);
	  print_results( \%res, $effects);
	  undef @vfs;# = ();	  
	  undef %res;# = ();
	  %grch37_remapping = ();
#	  exit if ($exit_count-- < 0);
	}
      }
    }
  }

  if ( @vfs || keys %res ) {
    my $effects = variation_effects(\@vfs) if ( ! $no_ensembl);
    print_results( \%res, $effects);
    @vfs = ();
    %res = ();
    %grch37_remapping = ();
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
    
    my @annotations = ('position', 'change', 'filter', 'score');
    push @annotations, ('depth', '', '', '', '', '') if ( $bam );
      
    push @annotations, ('gene', 'transcript', 'region', 'codon pos', 'AA change') if ( ! $no_ensembl);
    push @annotations, ('Grantham score') if ( ! $grantham);
    push @annotations, ('external ref') if ( ! $no_ensembl);
    push @annotations, ('dbsnp freq') if ( $use_local_dbsnp);
    push @annotations, ('dbsnp flags') if ( $use_local_dbsnp);
    push @annotations, ('HGMD') if ( $hgmd);


    push @annotations, 'pfam' if ( $pfam );
    push @annotations, 'regulation' if ( $regulation );

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
      push @line, $$mapping{$name}{filter};
      push @line, $$mapping{$name}{qual};
      
      if ( $$mapping{$name}{base_dist} ) {
	push @line, $$mapping{$name}{base_dist}{total};      
	push @line, $$mapping{$name}{base_dist}{type};      
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
      push @line, $$mapping{$name}{filter};
      push @line, $$mapping{$name}{qual};
      if ( $$mapping{$name}{base_dist} ) {
	push @line, $$mapping{$name}{base_dist}{total};      
	push @line, $$mapping{$name}{base_dist}{type};      	
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
	push @effect_line, $$snp_effect{ grantham } || "" if ( $grantham );

	$$snp_effect{ rs_number } = "<a href='$dbsnp_link$$snp_effect{ rs_number }'>$$snp_effect{ rs_number }</a>"
	    if ( $$snp_effect{ rs_number } && $html_out && $$snp_effect{ rs_number } =~ /rs\d+/ && $$snp_effect{ rs_number } !~ /href/);

	$$snp_effect{ rs_number } = "<a href='$ens_snp_link$$snp_effect{ rs_number }'>$$snp_effect{ rs_number }</a>"
	    if ( $$snp_effect{ rs_number } && $html_out && $$snp_effect{ rs_number } =~ /ENSSNP\d+/ && $$snp_effect{ rs_number } !~ /href/);

	push @effect_line, $$snp_effect{ rs_number   } || "";
	push @effect_line, $$snp_effect{ dbsnp_freq  } || "" if ( $use_local_dbsnp );
	push @effect_line, $$snp_effect{ dbsnp_flags } || "" if ( $use_local_dbsnp );

	push @effect_line, $$snp_effect{ HGMD } || "" if ( $hgmd );


#	die Dumper( $snp_effect ) if ( $$snp_effect{ dbsnp_freq } );

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
      push @line, $$mapping{$name}{base_dist}{type};      
      map { push @line, $$mapping{$name}{base_dist}{$_} if ($$mapping{$name}{base_dist}{$_})} ( 'A', 'C', 'G', 'T', 'N');
    }

    push @line, $$effect[0]{ rs_number } || "";
    push @line, $$effect[0]{ dbsnp_freq  } || "" if ( $use_local_dbsnp );
    push @line, $$effect[0]{ dbsnp_flags } || "" if ( $use_local_dbsnp );

    push @line, $$effect{ HGMD } || "" if ( $hgmd );

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

      push @effect_line, $$snp_effect{ rs_number };
      push @effect_line, $$snp_effect{ dbsnp_freq  } || "" if ( $use_local_dbsnp );
      push @effect_line, $$snp_effect{ dbsnp_flags } || "" if ( $use_local_dbsnp );

      push @effect_line, $$snp_effect{ HGMD } || "" if ( $hgmd );


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


#	$existing_vf = "";
    my ($dbsnp_flags, $dbsnp_freq) = ('','');
    if ( $use_local_dbsnp ) {
      
      my ($chr, $pos);
      if ($from_36 && $grch37_remapping{$vf->variation_name()}) {
	($chr, $pos) = split(":", $grch37_remapping{$vf->variation_name()});
#	print "looking for snp at: $chr $pos ($existing_vf) ". ($vf->variation_name())."\n";
      }
      else {
	($chr, $pos) = split(":", $vf->variation_name());
      }
     
      $sth_dbsnp->execute($chr, $pos);
      my $result = $sth_dbsnp->fetchrow_hashref();

      if ( $result->{rs} ) {
	$existing_vf = $result->{rs};
	$dbsnp_flags = $result->{flags} || "";
	$dbsnp_flags =~ s/RV;//;
	$dbsnp_flags =~ s/dbSNPBuildID=\d+;//;
	$dbsnp_flags =~ s/WGT=\d+;//;
	$dbsnp_flags =~ s/SLO;//;
	$dbsnp_flags =~ s/VC=\w+?;//;
	$dbsnp_flags =~ s/VP=.+?;//;
	
	
	$sth_pop->execute($existing_vf);
	
	my ( $observations, $seen, $freq, $pop) = (0,0, 0);
	while (my $ref2 = $sth_pop->fetchrow_hashref ) {
	  $pop++;
	  my $count = int($$ref2{allele_freq}*$$ref2{sample_size});
	  $observations += $$ref2{sample_size};
	  $seen += $count;
	  $freq += $$ref2{allele_freq};
	}
	$dbsnp_freq = sprintf("$seen/%.2f", $freq/$pop) if ( $pop);
#	print "($dbsnp_freq, $dbsnp_flags) =); \n";
	
      }
	  
    }

    my $HGMD;
    if ( $hgmd ) {
      $HGMD = "N";
      my ($chr, $pos) = split(":", $vf->variation_name());
      
#      my $sth = $hgmd->prepare("SELECT rs, chr, pos FROM snp WHERE (chr='$chr' AND pos = '$pos') OR rs='$existing_vf';");
      $sth_hgmd->execute($chr, $pos, $existing_vf);
      my $result = $sth_hgmd->fetchrow_hashref();
      
      $HGMD = "Y" if ( defined $result->{rs} || defined $result->{chr} );
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

	$gene_res{ regulation } = $regulations; 

	$gene_res{ rs_number }   = $existing_vf || "";
	$gene_res{ dbsnp_freq }  = $dbsnp_freq || "";
	$gene_res{ dbsnp_flags } = $dbsnp_flags || "";
	$gene_res{ HGMD } = $HGMD || "";

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
  my $largest = 0;
  foreach my $base (sort {$base_stats{$b} <=> $base_stats{$a}} keys %base_stats ) {
    my $perc = sprintf("%.2f", $base_stats{$base}/$total*100);
    $largest = $perc if ( $largest < $perc );
    my $qual = $qual_stats{$base} || 0;
    $res{$base} = "$base: $base_stats{$base}($perc%)/$qual";
  }

  $res{type}  = "homozygous"   if ( $largest > 76);
  $res{type}  = "heterozygous" if ( $largest < 75 && $largest > 35);
  $res{type}  = "unknown"      if ( $largest < 35);
  $res{type} .= " ($largest)";

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
				filter    => $filter,
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
# Kim Brugger (01 Dec 2010)
sub grantham_score {
  my ($aa1, $aa2) = @_;

  $aa1 = uc($aa1);
  $aa2 = uc($aa2);
  
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
  
  print "USAGE: $0 -b[am file] -v[cf file] -T<ranform, use if mapped against hg18> -H<tml out put> -I<GV links> -p<fam domains> -r<egulation information>\n";
  print "USAGE: -E[nsembl.ensembl.org databases, SLOOOOOOWWWWWERRRRRR ]\n";  
  print "USAGE: -g[rantham score]\n";  
  print "USAGE: -L[local SNP database backup, as ensembl is broken!!]\n";  
  print "USAGE: -m[in mapping qual. for base dist. report]\n";  
  print "USAGE: -q[uality cutoff for SNPs to report]\n";
  print "USAGE: -f[ multi line report]\n";
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
