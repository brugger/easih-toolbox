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

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 1;
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
my $leeway       = $opts{l} || 100;

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

if ( $bam && $baits ) {

  $DACstring = DepthAndCoverage($bam,$baits) if ( $baits);
  ($ANstring, my @SMstrings) = BamHead_ANSM($bam);

  if($ANstring eq "hg18" && ! $from_36) {
    $from_36 = 1;
    print STDERR "**************************************\n";
    print STDERR "\nSetting -T flag automatically as reference is hg18\nPlease use -T<ransform> option so that the coordinates can be found in Ensembl\n";
    print STDERR "**************************************\n";
  }

  die "\n\n ERROR!!!! Wrong bed file used: $baits !~ $ANstring \n\n", if($baits !~ m/$ANstring/);
  
  $ANstring = "# Aligned2Reference: $ANstring";
  
  $SMstring = "# SampleName: ";
  
  foreach my $string(@SMstrings) {
    $SMstring .= "$string,";
  }
  
  $SMstring =~ s/,$//;
}



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
readin_vcf( $snp_indel_vcf ) if ( $snp_indel_vcf);

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
    $res{$position}{base_dist} = base_dist( $chr, $pos, $res{$position}{ref_base}, $res{$position}{alt_base}) if ( $bam && $basecount );
    
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


my $printed_header = 0;

# 
# 
# 
# Kim Brugger (08 Jul 2010)
sub print_results {
  my ( $mapping ) = @_;

  my (@res, @filtered_res);
  if ( ! $printed_header ) {

    $0 =~ s/.*\///;
    my @header;
    push @header, [ "#EASIH Variation Report v1.20.1"];
    push @header, [ "#commandline: $0 @argv"];
    push @header, [ "#dbases: ". EASIH::SNPs::db_info()];
    push @header, [ "#script version: $version"];
    
    push @header, [ "#bait filtering with a leeway of: $leeway and $baits as the bait file"] if ($baits );
    
    push @header, ["$DACstring"] if ( $DACstring );
    push @header, ["$SMstring"]  if ( $SMstring  );
    push @header, ["$ANstring"]  if ( $ANstring  );
  
    push @header, [ q{##INFO=<Field=Position,ExampleValue=X:2715425,Description="Chromosome and coordinate">}];
    push @header, [ q{##INFO=<Field=Change,ExampleValue=A>G,Description="Base change">}];
    push @header, [ q{##INFO=<Field=Filter,ExampleValue=PASS,Description="SNPs can be filtered based on certain criteria such as sequence depth">}];
    push @header, [ q{##INFO=<Field=Score,ExampleValue=1548.43,Description="SNP confidence score determined by GATK">}];
    push @header, [ q{##INFO=<Field=Depth,ExampleValue=11,Description="Number of reads covering base">}];
    push @header, [ q{##INFO=<Field=Genotype,ExampleValue=HOMO,Description="Nature of SNP">}];
    push @header, [ q{##INFO=<Field=Gene,ExampleValue=XG,Description="Gene Name or Ensembl gene ID, if there is no gene name ">}];
    push @header, [ q{##INFO=<Field=Transcript,ExampleValue=NM_175569.2,Description="Refseq id or Ensembl transcript ID if no efseq id">}];
    push @header, [ q{##INFO=<Field=Effect,ExampleValue=NON_SYNONYMOUS_CODING,Description="Indicates whether the variation is in gene, etc">}];
    push @header, [ q{##INFO=<Field=Codon pos,ExampleValue=c.546,Description="Which codon the variation occurs in">}];
    push @header, [ q{##INFO=<Field=AA change,ExampleValue=p.Ser157 Phe,Description="Amino acid change: p.FromPosition To">}];
    push @header, [ q{##INFO=<Field=Grantham score,ExampleValue=155,Description="Measure of difference between reference and mutation amino acid">}];
    push @header, [ q{##INFO=<Field=dbsnp,ExampleValue=rs111382948,Description="External Database ref e.g., dbSNP">}];
    push @header, [ q{##INFO=<Field=dbsnp flags,ExampleValue=VLD;KGPilot1VC=SNP,Description="dbSNP details">}];
    push @header, [ q{##INFO=<Field=HGMD,ExampleValue=Y,Description="Is the SNP present in the HGMD database">"}];
    push @header, [ q{##INFO=<Field=pfam,ExampleValue=Sulfatase,Description="PFAM domain">}];
    push @header, [ q{##INFO=<Field=PolyPhen,ExampleValue=tolerated/0.59,Description="Estimates if the SNP is pathogenic">}];
    push @header, [ q{##INFO=<Field=SIFT,ExampleValue=benign/0,Description="Estimates if the SNP is pathogenic">}];
    push @header, [ q{##INFO=<Field=Condel,ExampleValue=neutral/0.391,Description="Estimates if the SNP is pathogenic">}];
    push @header, [ q{##INFO=<Field=GERP,ExampleValue=800.7,Description="The genomic evolutionary rate for the nucleotide in mammals">}];
    
    
    push @res, @header;
    push @filtered_res, @header;
    push @filtered_res, ["#filtered for most important snp effect"];
    
    my @annotations = ('Position', 'Change', 'Filter', 'Score', 'Depth', 'Genotype');
    push @annotations, ('', '', '', '', '') if ( $bam && $basecount);
    
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
    
    $printed_header++;  
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

    $$mapping{$name}{res} = sort_effects( $$mapping{$name}{res});

    my $first = 1;

    foreach my $effect ( @{$$mapping{$name}{res}} ) {

    # find the most important/interesting variation effect
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

    if (! $vf->slice->sub_Slice($vf->start, $vf->end, $vf->strand) ) {
      print STDERR "should this be hg18? cannot get slice: ". $vf->start .":". $vf->end ."-". $vf->strand . "\n";
      next;
    }
    my $cons = $ce_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($ce_mlss,$vf->slice->sub_Slice($vf->start, $vf->end, $vf->strand));

    # get consequences
    # results are stored attached to reference VF objects
    # so no need to capture return value here

    my $all_trans_var = $vf->get_all_TranscriptVariations;

    #print "LENGTH: [", scalar @$all_trans_var, "]\n";
    if (scalar @$all_trans_var == 0) {
	    push @{$$vars{ $name }{res}}, {};
    }

    foreach my $tv (@{$all_trans_var}) {

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
	      
 	      next if ( $logic_name ne 'pfam');
	      
 	      if ($tv->translation_start >= $prot_feat->start() and
 		  $tv->translation_end <= $prot_feat->end() ) {
		
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

    next if ((/^chromosome/) && (! $leeway)); 
    
    chomp;
    my ($chr, $start, $end) = split("\t", $_);

    ($chr, $start, $end) = $_ =~ /(.*?):(\d+)-(\d+)/
	if ( ! $start );

    next if ( ! $chr );

    $chr =~ s/chr//;
    
    if ($leeway)
    {
	$start -= $leeway;
	$end   += $leeway;
    }
    
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
# Kim Brugger (21 Jul 2010)
sub names_n_lengths {
  my ( $bam_file ) = @_;
  
  my $samtools  = find_program('samtools');
  my @sequences = ();
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);

    my ($name, $length);
    foreach my $field ( split("\t") ) {
      $name   = $1 if ( $field =~ /SN:(.*)/);
      $length = $1 if ( $field =~ /LN:(\d+)/);
    }

    push @sequences, [$name, $length] if ( $name && $length );
  }

    
  return \@sequences;
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


sub DepthAndCoverage
{
    my($bamfile,$baitfile) = @_;

    my $regions_ref  = readin_bed( $baitfile, 0 ) if ( $baitfile ); 
    my $samtools  = find_program('samtools');
    my $bam2depth = find_program('bam2depth');

    my $seqs = names_n_lengths( $bamfile );

    my ($START, $END) = (0, 1);
    
    my $total_reads;
    my %base_coverage;
    my %exon_coverage;
    
    my $escape = 1;
    
    my $patched_start = 0;
    
    foreach my $chr ( sort keys %$regions_ref ) {
	
#  $chr ="chr2";
	
	my (undef, $base_pos, $depth);
	my @regions =  @{$$regions_ref{$chr}};
	my $region  = shift @regions;
	my ($start, $end) = @{$region};
	
#  print "Doing $chr ... expecting ".@regions." regions\n";
	
	$chr = "chr$chr" if ($baits =~ /hg18/);

	open (my $pipe, "$bam2depth $bamfile $chr | " ) || die "Could not open pipe: $!\n";
	my $covered = 0;
	while(<$pipe>) {
	    chomp;
	    (undef, $base_pos, $depth) = split("\t", $_);
	    
	    next if ( $base_pos < $start);
	    
	    if ( $base_pos >= $start && 
		 $base_pos <= $end ) {
		
		$base_coverage{ $depth }++;
		$covered++;
	    }
	    else {
		
		if ( $start != $end ) {
		    
		    $exon_coverage{ int( $covered/($end-$start + 1)*100)}++;
		    $base_coverage{ 0 } += $end - $start + 1 - $covered;
		}
		
		($start, $end) = (undef, undef);
		$region  = shift @regions;
		last if ( ! $region );
		($start, $end) = @{$region};
		
		$covered = 0;
	    }      
	}
	close( $pipe );
	
	if (  $start ) {
	    if ( $start != $end ) {
		$exon_coverage{ int( $covered/($end-$start + 1)*100)}++;
		$base_coverage{ 0 } += $end - $start + 1 - $covered;
	    }
	}
	
	if (@regions) {
	    foreach $region (@regions) {
		($start, $end) = @{$region};
		$exon_coverage{ 0 }++;
		$base_coverage{ 0 } += $end - $start + 1;
	    }
	}
	
#  last;
    }
    
    my $DACstring;
    
    $DACstring = Depth($bamfile,%base_coverage);
    
    $DACstring .= Coverage(%exon_coverage);
    
      
    return($DACstring);
}


sub Depth
{
    my($bamfile,%base_coverage) = @_;
    my %binned = ( 0 => 0,
		   1 => 0,
		   10 => 0,
		   20 => 0,
		   30 => 0,
		   40 => 0);
    my $base_total = 0;
    foreach my $depth ( keys %base_coverage ) {  
	
	$base_total += $base_coverage{ $depth };
	if ( $depth == 0 || $depth == 1 ) {
	    $binned{ $depth } = $base_coverage{ $depth };
	    next;
	}
	elsif ( $depth < 10 ) {
	    $binned{ 10 } +=  $base_coverage{ $depth };
	}
	elsif ( $depth < 20 ) {
	    $binned{ 20 } +=  $base_coverage{ $depth };
	}
	elsif ( $depth < 30 ) {
	    $binned{ 30 } +=  $base_coverage{ $depth };
	}
	else {
	    $binned{ 40 } +=  $base_coverage{ $depth };
	}
    }

    my $Dstring;

    #$Dstring .= "# $bamfile\n";
    $Dstring .= "# base depth distribution:\n"; 
    $Dstring .=  sprintf("#     0x depth: $binned{  0 } (%.2f %%) \n", $binned{  0 }/$base_total*100);
    $Dstring .=  sprintf("#     1x depth: $binned{  1 } (%.2f %%) \n", $binned{  1 }/$base_total*100);
    $Dstring .= sprintf("#  2-10x depth: $binned{ 10 } (%.2f %%) \n", $binned{ 10 }/$base_total*100);
    $Dstring .= sprintf("# 10-20x depth: $binned{ 20 } (%.2f %%) \n", $binned{ 20 }/$base_total*100);
    $Dstring .= sprintf("# 20-30x depth: $binned{ 30 } (%.2f %%) \n", $binned{ 30 }/$base_total*100);
    $Dstring .= sprintf("#   >30x depth: $binned{ 40 } (%.2f %%) \n", $binned{ 40 }/$base_total*100);
    
    
    return($Dstring);
}


sub Coverage
{
    my %exon_coverage = @_;

    my %exons = (  0 => 0,
		   70 => 0,
		   80 => 0,
		   90 => 0,
		   100 => 0,
		   101 => 0);
    
    my $exon_total = 0;
    foreach my $coverage ( keys %exon_coverage ) {  
	
	$exon_total += $exon_coverage{ $coverage };
	if ( $coverage == 0 ) {
	    $exons{ $coverage } += $exon_coverage{ $coverage };
	}
	elsif ( $coverage < 70 ) {
	    $exons{ 70 } += $exon_coverage{ $coverage };
	}
	elsif ( $coverage < 80 ) {
	    $exons{ 80 } += $exon_coverage{ $coverage };
	}
	elsif ( $coverage < 90 ) {
	    $exons{ 90 } += $exon_coverage{ $coverage };
	}
	elsif ( $coverage < 100 ) {
	    $exons{ 100 } += $exon_coverage{ $coverage };
	}
	elsif ( $coverage >= 100 ) {
	    $exons{ 101 } += $exon_coverage{ $coverage };
	}
    }

    my $Cstring;

    $Cstring .= "# bait regions coverage:\n"; 
    $Cstring .= sprintf("#      0 %% of bait covered: $exons{  0 } (%.2f %%) \n", $exons{     0 }/$exon_total*100);
    $Cstring .= sprintf("#    <70 %% of bait covered: $exons{  70 } (%.2f %%) \n", $exons{   70 }/$exon_total*100);
    $Cstring .= sprintf("#  70-80 %% of bait covered: $exons{  80 } (%.2f %%) \n", $exons{   80 }/$exon_total*100);
    $Cstring .= sprintf("#  80-90 %% of bait covered: $exons{  90 } (%.2f %%) \n", $exons{   90 }/$exon_total*100);
    $Cstring .= sprintf("#  90-99 %% of bait covered: $exons{  100 } (%.2f %%) \n", $exons{  100 }/$exon_total*100);
    $Cstring .= sprintf("#    100 %% of bait covered: $exons{  101 } (%.2f %%)", $exons{  101 }/$exon_total*100);

    return($Cstring);
}


sub BamHead_ANSM 
{
  my ($bam_file) = @_;
  
  my $samtools  = find_program('samtools');
 
  my($ANstring,$SMstring) = ("Undef","Undef");  
  my(@SMstrings);

  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  
  my $count = 0;
  
  while(<$spipe>) 
  {
      if(/\@RG/)
      {
	  if(/SM:(.*?)\s/)
	  {
	      $SMstring = "$1";
	  }
	  
	  #$SMstring = "# SampleName: $SMstring";
	  push @SMstrings, $SMstring;
      }
            
      next, if($count);
      
      if(/\@SQ/)
      {
	  #$ANstring = "$1", if(/AN:(.*?)\s*/);
	  $ANstring = "$1", if(/AN:(.*)\s\z/);
	  $count++;
      }
  }
  
  #$ANstring = "# Aligned2Reference: $ANstring";

  $ANstring = "hg19" if ( $ANstring eq "human_g1k_v37");

  return($ANstring,@SMstrings);
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



