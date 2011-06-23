package EASIH::QC::db;
#
# 
# 
# 
# Kim Brugger (23 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use DBI;
my $dbi;




# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_qvs {
  my ($fid) = @_;
  
  my $q = "SELECT x,q0,q1,q2,q3,q4 FROM qv_boxplot where fid = ? ORDER BY x";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid ) || die "$DBI::errstr";
  my @res;
  while ( my @line = $sth->fetchrow_array()) {
    push @res, [@line];
  }

  return @res;
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_qvs {
  my ( $fid, $values ) = @_;

  foreach my $value ( @$values ) {
    _add_qv_set( $fid, @$value);
  }
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub _add_qv_set {
  my ( $fid, $x, $q0, $q1, $q2, $q3, $q4 ) = @_;

  my $q = "INSERT INTO qv_boxplot (fid, x, q0, q1, q2, q3, q4) VALUES (?,?,?,?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid, $x, $q0, $q1, $q2, $q3, $q4 ) || die "$DBI::errstr";
}




# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_file {
  my ( $file, $sample, $project, $platform, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC ) = @_;
  
  my $file_id = fetch_file_id( $file );
  return $file_id if ( $file_id );

  # check to see if it already exists.
  my $sample_id = fetch_sample_id( $sample);
  if ( !$sample_id ) {
    $sample_id = add_sample( $sample, $project);
  }

  my $q = "INSERT INTO file (sid, name, platform, sample_size, Q30bases, duplicates, partial_adaptors, Avg_AC) VALUES (?,?,?,?,?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $sample_id, $file, $platform, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC ) || die "$DBI::errstr";

  return $dbi->last_insert_id(undef, undef, qw(file fid));
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_file_id {
  my ( $file ) = @_;
  print "$file\n";
  my $q = "SELECT fid FROM file where name = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $file ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return $line[0] || undef;
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_sample {
  my ( $sample, $project ) = @_;

  # check to see if it already exists.
  my $sample_id = fetch_sample_id( $sample);
  return $sample_id if ($sample_id);

  
  # check to see if it already exists.
  my $project_id = fetch_project_id( $project);
  if ( !$project_id ) {
    $project_id = add_project( $project);
  }

  my $q = "INSERT INTO sample (pid, name) VALUES (?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $project_id, $sample ) || die "$DBI::errstr";

  return $dbi->last_insert_id(undef, undef, qw(sample sid));
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_sample_id {
  my ( $sample ) = @_;
  my $q = "SELECT sid FROM sample where name = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $sample ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return $line[0] || undef;
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_project {
  my ( $project ) = @_;
  
  # check to see if it already exists.
  my $project_id = fetch_project_id( $project);
  return $project_id if ( $project_id );

  my $q = "INSERT INTO project (name) VALUES (?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $project ) || die "$DBI::errstr";

  return $dbi->last_insert_id(undef, undef, qw(project pid));
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_project_id {
  my ( $project ) = @_;
  my $q = "SELECT pid FROM project where name = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $project ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return $line[0] || undef;
}



# 
# 
# 
# Kim Brugger (12 May 2011)
sub Connect {
  my (undef, $dbname, $dbhost) = @_;
  $dbname ||= "qc_dev";
  $dbhost ||= "mgpc17";

  $dbi = DBI->connect("DBI:mysql:$dbname:$dbhost", 'easih_admin', 'easih') || die "Could not connect to database: $DBI::errstr";


}


BEGIN {
  Connect();
}


1;



__END__







