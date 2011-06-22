package EASIH::SNPs;
#
# Interface to local dbsnp + HGMD database + phylop_phast
# 
# 
# 
# Kim Brugger (08 Dec 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;

use DBI;

my $snp_db;
my $pp_db;
my $dbi;
my $ppdbi;
my $ref_id_hg18;
my $ref_id_GRCh37;

my $sth_fetch_snp = $dbi->prepare("SELECT snp.*, ref.name as ref_name FROM snp, reference as ref WHERE chr=? AND pos = ? AND snp.ref_id=? AND snp.ref_id=ref.ref_id");
my $sth_fetch_rs  = $dbi->prepare("SELECT snp.*, ref.name as ref_name FROM snp, reference as ref WHERE rs=? AND snp.ref_id=ref.ref_id");
my $sth_fetch_pop = $dbi->prepare("select * from population where rs=?;");

my %populations;

my $sth_phylop = $ppdbi->prepare("SELECT score from phylop where chr=? AND pos = ? AND ref_id=?");
my $sth_phast  = $ppdbi->prepare("SELECT score from phast  where chr=? AND pos = ? AND ref_id=?");


# 
# 
# 
# Kim Brugger (12 Jan 2011)
sub db_info {

  return "connected to $snp_db at mgpc17";
  
}


#
#
#
# Kim Brugger (08 Dec 2010)
sub phast_score_GRCh37 {
  my ( $chr, $pos) = @_;
  return undef if ( ! $chr || !$pos);

  $sth_phast->execute( $chr, $pos, $ref_id_GRCh37);
  my $result = $sth_phast->fetchrow_hashref();

  return $result->{score} || '';
}


#
#
#
# Kim Brugger (08 Dec 2010)
sub phast_score_hg18 {
  my ( $chr, $pos) = @_;
  return undef if ( ! $chr || !$pos);

  $sth_phast->execute( $chr, $pos, $ref_id_hg18);
  my $result = $sth_phast->fetchrow_hashref();

  return $result->{score} || '';
}



# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub phylop_score_GRCh37 {
  my ( $chr, $pos) = @_;
  return undef if ( ! $chr || !$pos);
  
  $sth_phylop->execute( $chr, $pos, $ref_id_GRCh37);
  my $result = $sth_phylop->fetchrow_hashref();
  
  return $result->{score} || '';
}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub phylop_score_hg18 {
  my ( $chr, $pos) = @_;
  return undef if ( ! $chr || !$pos);
  
  $sth_phylop->execute( $chr, $pos, $ref_id_hg18);
  my $result = $sth_phylop->fetchrow_hashref();
  
  return $result->{score} || '';
}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub fetch_populations {
  
  my $sth = $dbi->prepare("select * from populations;");
  $sth->execute();
  while (my $res = $sth->fetchrow_hashref()) {
    $populations{ $res->{pop_id} }{short} = $res->{short};
    $populations{ $res->{pop_id} }{descr} = $res->{descr};
  }

}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub population_stats {
  my ( $rs, $population) = @_;

  $sth_fetch_pop->execute( $rs );
	
  my ( $observations, $seen, $freq, $pop) = (0,0, 0);
  while (my $ref2 = $sth_fetch_pop->fetchrow_hashref ) {
    $pop++;
    my $count = int($$ref2{allele_freq}*$$ref2{sample_size});
    $observations += $$ref2{sample_size};

    use Data::Dumper;

    return "$count/$$ref2{allele_freq}" if ( $population && $populations{ $ref2->{pop_id} }{short} eq $population);

    $seen += $count;
    $freq += $$ref2{allele_freq};
  }

  return "" if ( $population);

  return sprintf("$seen/%.2f", $freq/$pop) if ( $pop);
  return "";
}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub fetch_flags {
  my ( $rs) = @_;
  return undef if ( ! $rs);

  my @results;
  $sth_fetch_rs->execute( $rs);
  
  my $result = $sth_fetch_rs->fetchrow_hashref();

  return undef if ( ! $result->{flags} );
  my $flags = $result->{flags};
  $flags =~ s/RV;//;
#  $flags =~ s/dbSNPBuildID=\d+;//;
  $flags =~ s/WGT=\d+;//;
  $flags =~ s/SLO;//;
  $flags =~ s/VC=\w+?;//;
  $flags =~ s/VP=.+?;//;
  return $flags;
  
}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub fetch_rs {
  my ( $rs) = @_;
  return undef if ( ! $rs);

  my @results;
  $sth_fetch_rs->execute( $rs);
  
  while (my $result = $sth_fetch_rs->fetchrow_hashref() ) {
    push @results, $result;
  }
  
  return \@results;
}
 

# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub fetch_snp_hg18 {
  my ( $chr, $pos) = @_;
  return undef if ( ! $chr || !$pos);
  
  $sth_fetch_snp->execute( $chr, $pos, $ref_id_hg18);
  my $result = $sth_fetch_snp->fetchrow_hashref();
  
  return $result;
}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub fetch_snp_GRCh37 {
  my ( $chr, $pos) = @_;
  return undef if ( ! $chr || !$pos);
  
  $sth_fetch_snp->execute( $chr, $pos, $ref_id_GRCh37);
  
  my $result = $sth_fetch_snp->fetchrow_hashref();
#  while ( my $res2  = $sth_fetch_snp->fetchrow_hashref()) {
#    $$result{rs} .= ";". $$res2{rs};
#    $$result{hgmd} = $$res2{hgmd} if ($$res2{hgmd} eq "Y");
#  }

  return $result;
}

# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub fetch_ref_id {
  my ($name) = @_;

  my $sth = $dbi->prepare("select ref_id from reference where name='$name';");
  $sth->execute();
  my @entries =  $sth->fetchrow_array();
  if (@entries) {
    return $entries[0];
  }

  die "Could not find or insert reference '$name'\n";
}

# 
# Hard coded for now, change later...
# 
# Kim Brugger (08 Dec 2010)
BEGIN {
  $snp_db = "dbsnp_132";
  $pp_db  = "phylop_phast";

  $dbi = DBI->connect("DBI:mysql:$snp_db:mgpc17", 'easih_ro') || die "Could not connect to database: $DBI::errstr";
  $ppdbi = DBI->connect("DBI:mysql:$pp_db:mgpc17", 'easih_ro') || die "Could not connect to database: $DBI::errstr";
  $ref_id_hg18   = fetch_ref_id('hg18');
  $ref_id_GRCh37 = fetch_ref_id('GRCh37');
  fetch_populations();
}


1;
