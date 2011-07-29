package EASIH::SNPs;
#
# Interface to the local dbsnp
# 
# 
# 
# Kim Brugger (08 Dec 2010), contact: kim.brugger@easih.ac.uk
#
# Revised 12 May 2011, kim.brugger@easih.ac.uk
# 


use strict;
use warnings;

use DBI;

my $snp_db;
my $pp_db;
my $dbi;
my $ppdbi;

my $sth_fetch_snp;
my $sth_fetch_rs;

my $dbhost;
my $dbname;


# 
# 
# 
# Kim Brugger (12 May 2011)
sub New {
  (undef, $dbname, $dbhost) = @_;
  $dbname ||= "dbsnp_132_human";
  $dbhost ||= "mgpc17";

  $dbi = DBI->connect("DBI:mysql:$dbname:$dbhost", 'easih_ro') || die "Could not connect to database: $DBI::errstr";

  $sth_fetch_snp = $dbi->prepare("SELECT snp.* FROM snp WHERE chr=? AND pos = ?");
  $sth_fetch_rs  = $dbi->prepare("SELECT snp.* FROM snp WHERE rs=?");

}

# 
# 
# 
# Kim Brugger (12 Jan 2011)
sub db_info {
  return "connected to $dbname at $dbhost";
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
  return $result->{flags};
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
sub fetch_snp {
  my ( $chr, $pos) = @_;
  return undef if ( ! $chr || !$pos);
  
  $sth_fetch_snp->execute( $chr, $pos );
  
  my $result = $sth_fetch_snp->fetchrow_hashref();
  return $result;
}

<<<<<<< HEAD
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
