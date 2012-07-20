package EASIH::Config;
# 
# Global configuration stuff
# 
# 
# Kim Brugger (05 Jul 2011), contact: kim.brugger@easih.ac.uk


use vars qw ( @illumina_dirs
              %libraries
              $test_db_host
              $test_db_user
              $test_db_passwd);




$test_db_host    = "mgpc17.medschl.ac.uk";
$test_db_user    = "easih_admin";
$test_db_passwd  = "easih";



# 
# 
# 
# Kim Brugger (27 Feb 2012)
sub db_config {
  return ($test_db_host, $test_db_user, $test_db_passwd);
}

@illumina_dirs = ('/seqs/illumina2/',
		  '/seqs/babraham/',
		  '/seqs/illumina3',  
		  '/seqs/illumina4',  
		  '/seqs/illumina5');

%libraries = (
  Human          => "/data/refs/human_1kg/bowtie/human_g1k_v37",
  NCBI_ref       => "/data/refs/archive/human_genomic_transcript/ref_contig",
  NCBI_alt       => "/data/refs/archive/human_genomic_transcript/alt_contig_HuRef",
  NCBI_rna       => "/data/refs/archive/human_genomic_transcript/rna",
  Mouse          => "/data/refs/mm9/bowtie/mm9",
  Rat            => "/data/refs/archive/rn4/bowtie/rn4",
  Ecoli          => "/data/refs/archive/Ecoli/U00096_2",
  MyRDB_bact     => "/data/refs/archive/bacteria/myRDP-bacteria",
  GreenGenes     => "/data/refs/archive/greengenes/greengenes_unaligned",
  Yeast          => "/data/refs/archive/Scerevisiae/yeast",
  Pombe          => "/data/refs/archive/pombe/pombe",
  Viruses        => "/data/refs/archive/viruses/gb_virus_seq",
  PhiX           => "/data/refs/archive/PhiX/PhiX",
  Adapters       => "/data/refs/archive/fastqc_contaminants/Contaminants",
  Vectors        => "/data/refs/archive/UniVec/UniVec",
  Pichia         => "/data/refs/archive/pichia/SO/bowtie/pichia",
  Trypanosoma    => "/data/refs/archive/trypanosoma/Tb927_genome_230210",
  Tryp_gamb      => "/data/refs/archive/trypanosoma/Tbgamb_02_v2",
  Leishmania     => "/data/refs/archive/leishmania/leishmania",
  Ribo_ARB_SSU   => "/data/refs/archive/ribo/hs_ssu_r106_embl",
  Ribo_ARB_LSU   => "/data/refs/archive/ribo/hs_lsu_r106_embl",
  Ribo_UNI       => "/data/refs/archive/ribo/mart_uniprot_ribosomal",
  Ribo_ENS       => "/data/refs/archive/ribo/ens_all_ribo",
  Ribo_GO        => "/data/refs/archive/ribo/amigo_ens_anno_mf",
  MART_RNA_GENES => "/data/refs/archive/ribo/mart_rna_genes",
  MART_RRNA_ONLY => "/data/refs/archive/ribo/mart_rrna_only",
  MART_RRNA_PLUS => "/data/refs/archive/ribo/mart_rrna_mtrrna_pseudorrna",);



1;
