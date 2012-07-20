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
use Data::Dumper;

use DBI;

my $snp_db;
my $pp_db;
my $dbi;
my $ppdbi;

my $sth_fetch_snp;
my $sth_fetch_rs;
my $sth_fetch_cm_up;
my $sth_fetch_cm_down;

my $dbhost;
my $dbname;


# 
# 
# 
# Kim Brugger (12 May 2011)
sub New {
  (undef, $dbname, $dbhost) = @_;
  $dbname ||= "dbsnp_135";
  $dbhost ||= "mgpc17.medschl.cam.ac.uk";

  $dbi = DBI->connect("DBI:mysql:$dbname:$dbhost", 'easih_ro') || die "Could not connect to database: $DBI::errstr";

  $sth_fetch_snp      = $dbi->prepare("SELECT snp.* FROM snp WHERE chr=? AND pos = ?");
  $sth_fetch_rs       = $dbi->prepare("SELECT snp.* FROM snp WHERE rs=?");
  $sth_fetch_cm_up    = $dbi->prepare("SELECT snp.* FROM snp WHERE chr = ? and pos <= ? AND centimorgan order by pos DESC limit 10; ");
  $sth_fetch_cm_down  = $dbi->prepare("SELECT snp.* FROM snp WHERE chr = ? and pos >= ? AND centimorgan order by pos limit 10; ");

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



# 
# 
# 
# Kim Brugger (16 Nov 2011)
sub CM {
  my ( $chr, $pos) = @_;

  return undef if ( ! $chr || !$pos);

  my @cms;

  $sth_fetch_cm_up->execute( $chr, $pos );
  while (my $result = $sth_fetch_cm_up->fetchrow_hashref() ) {
    unshift @cms, [$$result{pos}, $$result{centimorgan}];
  }

  $sth_fetch_cm_down->execute( $chr, $pos );
  while (my $result = $sth_fetch_cm_down->fetchrow_hashref() ) {
    push @cms, [$$result{pos}, $$result{centimorgan}];
  }
#  print "chr:$pos\n";
#  print Dumper( \@cms );
  
  for(my $i = 0; $i<@cms-1; $i++ ){ 
    return $cms[$i][1] if ($cms[$i][0] == $pos);
    if ($cms[$i][0] < $pos && $cms[$i + 1][0] > $pos) {
      return $cms[$i + 1][1] if ( $cms[$i][0] ==  $cms[$i + 1][0]);
      
      return sprintf("%.2f", ($cms[$i][1]+$cms[$i + 1][1])/2);
    }
  }

  
  return "NA";
}


1;
