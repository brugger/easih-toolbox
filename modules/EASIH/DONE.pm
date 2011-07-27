package EASIH::DONE;
# 
# Data-Offloading-aNd-Evaluation...
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Time::HiRes;

use DBI;

my $dbi;

# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub Connect {
  my ($dbname, $dbhost) = @_;
  $dbname ||= "done";
  $dbhost ||= "mgpc17.medschl.cam.ac.uk";

  $dbi = DBI->connect("DBI:mysql:$dbname:$dbhost", 'easih_admin','easih') || die "Could not connect to database: $DBI::errstr";
}



BEGIN {
  Connect();
}

# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_files_from_run {
  my ($run) = @_;
  my $sth = $dbi->prepare("SELECT f.name FROM file as f, run as r WHERE r.name = ? and r.rid = f.rid");
  $sth->execute( $run );
  my @res;
  while ( my @line =  $sth->fetchrow_array()) {
    push @res, @line;
  }

  return @res;
}


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_samples_from_run {
  my ($run) = @_;
  my $sth = $dbi->prepare("SELECT f.sid FROM file as f, run as r WHERE r.name = ? and r.rid = f.rid  group by f.sid;");
  $sth->execute( $run );
  my @res;
  while ( my @line =  $sth->fetchrow_array()) {
    push @res, @line;
  }

  return @res;
}


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_files_from_rid {
  my ($rid) = @_;
  my $sth = $dbi->prepare("SELECT name FROM file WHERE rid = ?");
  $sth->execute( $rid );
  my @res;
  while ( my @line =  $sth->fetchrow_array()) {
    push @res, @line;
  }

  return @res;
}


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_files_from_rid_by_sample {
  my ($rid) = @_;
  my $sth = $dbi->prepare("select f.name, s.name from file as f, run as r, sample as s WHERE r.rid=? and r.rid = f.rid AND s.sid=f.sid");
  $sth->execute( $rid );
  my %res;
  while ( my @line =  $sth->fetchrow_array()) {
    push @{$res{ $line[1] }}, $line[0];
  }

  return %res;
}

# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_files {
  my $sth = $dbi->prepare("SELECT * FROM file");
  $sth->execute(  );
  my @res;
  while ( my @line =  $sth->fetchrow_array()) {
    push @res, \@line;
  }

  return @res;
}


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub run_log {
  my ($run ) = @_;

  my $rid = fetch_run_id( $run );

  my $entries  = fetch_offloading_entries($rid);
  my @files    = fetch_files_from_rid($rid);

  my $printstring;
  
  $printstring .= "========================================\n";
  $printstring .= "  $run\n";
  $printstring .= "========================================\n";



  use POSIX qw( strftime );
  foreach my $entry ( @$entries ) {
    my $microsec = $$entry[0] % 100000;
    $$entry[0] /= 100000;
    $$entry[0] = strftime("%d/%m/%Y %H:%M:%S.$microsec", localtime( $$entry[0] ));

    $printstring .=   "$$entry[0]\t$$entry[1]\n";
  }


  $printstring .=   "========================================\n";
  foreach my $file ( @files ) {

      $printstring .=    "$file\n";
  }

  $printstring .=   "========================================\n";

  return($printstring);
}



# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub add_offloading_status {
  my ($rid, $status) = @_;
  my $timestamp = Time::HiRes::gettimeofday()*100000;
  my $sth = $dbi->prepare("INSERT INTO offloading (rid, stamp, status) VALUES (?,?,?) ");
  $sth->execute( $rid, $timestamp, $status );
}

# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_offloading_entries {
  my ($rid) = @_;
  
  my $q = "SELECT * FROM offloading where rid = ? ORDER BY stamp";
  my $sth = $dbi->prepare($q);
  $sth->execute( $rid );
  
  my @res;
  while ( my $results = $sth->fetchrow_hashref()) {
    push @res, [$$results{stamp}, $$results{status}];
  }

  @res = sort { $$a[0] <=> $$b[0] } sort @res;

  return \@res;
}



# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_latest_offloading_status {
  my ( $rid ) = @_;

  my $entries = fetch_offloading_entries( $rid );
  my $status;
  $status =  $$entries[-1][1] if ($$entries[-1]);
  return $status;
}


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_platform {
  my ( $run ) = @_;

  my $q = "SELECT platform from run where name = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $run );

  my $platform = undef;
  my $results = $sth->fetchrow_arrayref();
  $platform = $$results[0] if ($results);
  return $platform;
}

  
# 
# 
# 
# Kim Brugger (22 Jun 2011)
sub fetch_offloading_failure {
  my ($rid) = @_;

  my $entries = fetch_offloading_entries( $rid );
  my $status;
  $status =  $$entries[-2][1] if ($$entries[-2]);
  return $status;
}

# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_adaptors {
  my ($fid) = @_;
  
  my $q = "SELECT x, percent FROM adaptors where fid = ? ORDER BY x";
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
sub add_adaptors {
  my ( $fid, $values ) = @_;

  my $check_file_name = fetch_filename( $fid );
  if ( ! $check_file_name ) {
    print "'$fid' does not exist in the database\n";
    return 0;
  }

  foreach my $value ( @$values ) {
    _add_adaptor( $fid, @$value);
  }
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub _add_adaptor {
  my ( $fid, $x, $perc ) = @_;


  my $q = "REPLACE INTO adaptors (fid, x, percent) VALUES (?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid, $x, $perc ) || die "$DBI::errstr";
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_duplicate_seqs {
  my ($fid) = @_;
  
  my $q = "SELECT sequence, percentage, source FROM duplicated_seqs where fid = ? ORDER BY percentage";
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
sub add_duplicate_seqs {
  my ( $fid, $values ) = @_;

  my $check_file_name = fetch_filename( $fid );
  if ( ! $check_file_name ) {
    print "'$fid' does not exist in the database\n";
    return 0;
  }

  foreach my $value ( @$values ) {
    _add_duplicate_seq( $fid, @$value);
  }
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub _add_duplicate_seq {
  my ( $fid, $sequence, $perc, $source ) = @_;

  my $q = "REPLACE INTO duplicated_seqs (fid, sequence, percentage, source) VALUES (?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute(  $fid, $sequence, $perc, $source ) || die "$DBI::errstr";
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_duplicates {
  my ($fid) = @_;
  
  my $q = "SELECT x,observations FROM duplicates where fid = ? ORDER BY x";
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
sub add_duplicates {
  my ( $fid, $values ) = @_;

  my $check_file_name = fetch_filename( $fid );
  if ( ! $check_file_name ) {
    print "'$fid' does not exist in the database\n";
    return 0;
  }

  foreach my $value ( @$values ) {
    _add_duplicate( $fid, @$value);
  }
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub _add_duplicate {
  my ( $fid, $x, $dups ) = @_;

  my $q = "REPLACE INTO duplicates (fid, x, observations) VALUES (?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid, $x, $dups ) || die "$DBI::errstr";
}

# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_gc_distribution {
  my ($fid) = @_;
  
  my $q = "SELECT x,percent_gc FROM GC_distribution where fid = ? ORDER BY x";
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
sub add_gc_distribution {
  my ( $fid, $values ) = @_;

  my $check_file_name = fetch_filename( $fid );
  if ( ! $check_file_name ) {
    print "'$fid' does not exist in the database\n";
    return 0;
  }

  foreach my $value ( @$values ) {
    _add_gc_dist( $fid, @$value);
  }
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub _add_gc_dist {
  my ( $fid, $x, $perc ) = @_;

  my $q = "REPLACE INTO GC_distribution (fid, x, percent_gc) VALUES (?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid, $x, $perc ) || die "$DBI::errstr";
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_base_distribution {
  my ($fid) = @_;
  
  my $q = "SELECT x,A,C,G,T,N FROM base_distribution where fid = ? ORDER BY x";
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
sub add_basedists {
  my ( $fid, $values ) = @_;

  my $check_file_name = fetch_filename( $fid );
  if ( ! $check_file_name ) {
    print "'$fid' does not exist in the database\n";
    return 0;
  }

  foreach my $value ( @$values ) {
    _add_basedist( $fid, @$value);
  }
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub _add_basedist {
  my ( $fid, $x, $A, $C, $G, $T, $N ) = @_;

  my $q = "REPLACE INTO base_distribution (fid, x, A, C, G, T, N) VALUES (?,?,?,?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid, $x, $A, $C, $G, $T, $N ) || die "$DBI::errstr";
}






# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_qvs_histogram {
  my ($fid) = @_;
  
  my $q = "SELECT x,height FROM qv_histogram where fid = ? ORDER BY x";
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
sub add_qvs_histogram {
  my ( $fid, $values ) = @_;

  my $check_file_name = fetch_filename( $fid );
  if ( ! $check_file_name ) {
    print "'$fid' does not exist in the database\n";
    return 0;
  }

  foreach my $value ( @$values ) {
    _add_qv_hist( $fid, @$value);
  }
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub _add_qv_hist {
  my ( $fid, $x, $height ) = @_;

  my $q = "REPLACE INTO qv_histogram (fid, x, height) VALUES (?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid, $x, $height ) || die "$DBI::errstr";
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_qvs_boxplot {
  my ($fid) = @_;
  
  my $q = "SELECT x,q0,q1,q2,q3,q4,q5,q6,q7 FROM qv_boxplot where fid = ? ORDER BY x";
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
sub add_qvs_boxplot {
  my ( $fid, $values ) = @_;

  my $check_file_name = fetch_filename( $fid );
  if ( ! $check_file_name ) {
    print "'$fid' does not exist in the database\n";
    return 0;
  }

  foreach my $value ( @$values ) {
    _add_qv_boxplot( $fid, @$value);
  }
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub _add_qv_boxplot {
  my ( $fid, $x, $q0, $q1, $q2, $q3, $q4, $q5, $q6, $q7 ) = @_;

  my $q = "REPLACE INTO qv_boxplot (fid, x, q0, q1, q2, q3, q4, q5, q6, q7) VALUES (?,?,?,?,?,?,?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid, $x, $q0, $q1, $q2, $q3, $q4, $q5, $q6, $q7 ) || die "$DBI::errstr";
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_mapping_stats {
  my ( $fid1, $fid2, $reference, $unique_hits, $nonunique_hits, $duplicates ) = @_;

  
  my $q = "REPLACE INTO mapping_stats ( fid1, fid2, reference, unique_hits, non_unique_hits, duplicates) VALUES (?,?,?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute(  $fid1, $fid2, $reference, $unique_hits, $nonunique_hits, $duplicates ) || die "$DBI::errstr";
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_mapping_stats {
  my ( $fid1, $fid2 ) = @_;
  my $sth;
  if ( $fid2 ) {
    my $q = "SELECT * FROM  mapping_stats where fid1 = ? AND fid2 = ?  ORDER BY reference";
    $sth = $dbi->prepare($q);
    $sth->execute( $fid1, $fid2 ) || die "$DBI::errstr";
  }
  else {
    my $q = "SELECT reference, unique_hits, non_unique_hits, duplicates FROM  mapping_stats where fid1 = ? OR fid2 = ? ORDER BY reference";
    $sth = $dbi->prepare($q);
    $sth->execute( $fid1, $fid1 ) || die "$DBI::errstr";
  }
  my @res;
  while(my @line =  $sth->fetchrow_array()) {
    push @res, [@line];
  }
  return @res;
}




# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_illumina_multiplex_stats {
  my ( $rid, $fid, $lane, $sample, $bcode, $total_reads, $pass_filter, $ratio ) = @_;

  my $q = "REPLACE INTO illumina_multiplex_stats ( rid, fid, lane, sample, bcode, total_reads, pass_filter, ratio) VALUES (?,?,?,?,?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $rid, $fid, $lane, $sample, $bcode, $total_reads, $pass_filter, $ratio ) || die "$DBI::errstr";
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_illumina_multiplex_stats_by_rid {
  my ( $rid ) = @_;
  my $q = "SELECT * FROM  illumina_multiplex_stats where rid = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $rid ) || die "$DBI::errstr";
  my @res;
  while(my @line =  $sth->fetchrow_array()) {
    push @res, [@line];
  }

  return @res;
}




# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_illumina_lane_stats {
  my ( $rid, $fid, $lane, $read_nr, $sample, $total_reads, $pass_filter ) = @_;

  my $q = "REPLACE INTO illumina_lane_stats ( rid, fid, lane, read_nr, sample, total_reads, pass_filter) VALUES (?,?,?,?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $rid, $fid, $lane, $read_nr, $sample, $total_reads, $pass_filter ) || die "$DBI::errstr";
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_illumina_lane_stats_by_rid {
  my ( $rid ) = @_;
  my $q = "SELECT * FROM  illumina_lane_stats where rid = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $rid ) || die "$DBI::errstr";
  my @res;
  while(my @line =  $sth->fetchrow_array()) {
    push @res, [@line];
  }
  return @res;
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub update_file {
  my ( $fid, $total_reads, $read_length, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC ) = @_;
  
  return if (!$total_reads && !$read_length && !$sample_size && !$Q30bases && !$duplicates && !$partial_adaptors & !$Avg_AC);

  my $q = "UPDATE file SET ";
  my @updates;
  push @updates, "total_reads = '$total_reads' "           if ( $total_reads );
  push @updates, "read_length = '$read_length' "            if ( $read_length );
  push @updates, "sample_size = '$sample_size' "           if ( $sample_size );
  push @updates, "Q30bases = '$Q30bases' "                 if ( $Q30bases );
  push @updates, "duplicates = '$duplicates' "             if ( $duplicates );
  push @updates, "partial_adaptors = '$partial_adaptors' " if ( $partial_adaptors );
  push @updates, "Avg_AC = '$Avg_AC' "                     if ( $Avg_AC );

  $q .= join(", ", @updates) . " WHERE fid='$fid'";

  my $sth = $dbi->prepare($q);
  $sth->execute( ) || die "$DBI::errstr";

}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_file {
  my ( $file, $sample, $project, $run, $platform, $total_reads, $read_length, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC ) = @_;

  my $file_id = fetch_file_id( $file );
  return $file_id if ( $file_id );

  # check to see if it already exists.
  my $sample_id = fetch_sample_id( $sample);
  if ( !$sample_id ) {
    $sample_id = add_sample( $sample, $project);
  }

  my ($run_id, undef) = fetch_run_id( $run);
  if ( !$run_id ) {
    $run_id = add_run( $run, $platform);
  }

  my $q = "REPLACE INTO file (sid, rid, name, total_reads, read_length, sample_size, Q30bases, duplicates, partial_adaptors, Avg_AC) VALUES (?,?,?,?,?,?,?,?,?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $sample_id, $run_id, $file, $total_reads, $read_length, $sample_size, $Q30bases, $duplicates, $partial_adaptors, $Avg_AC ) || die "$DBI::errstr";

  return $dbi->last_insert_id(undef, undef, qw(file fid));
}




# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_file_info {
  my ( $fid ) = @_;
  my $q = "SELECT name, total_reads, read_length, sample_size, Q30bases, duplicates, partial_adaptors, Avg_AC, sid, rid  FROM file where fid = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return @line;
}
  


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_filename {
  my ( $fid ) = @_;
  my $q = "SELECT name FROM file where fid = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $fid ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return $line[0] || undef;
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_file_id {
  my ( $file ) = @_;
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

  my $q = "REPLACE INTO sample (pid, name) VALUES (?,?)";
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
sub fetch_sample_name {
  my ( $sid ) = @_;
  my $q = "SELECT name FROM sample where sid = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $sid ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return $line[0] || undef;
}

# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_samples {
  my $q = "SELECT * FROM sample";
  my $sth = $dbi->prepare($q);
  $sth->execute( ) || die "$DBI::errstr";
  my @data;
  while (my @line =  $sth->fetchrow_array()) {
    push @data, [@line];
  }

  return @data;
}


# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub add_run {
  my ( $run, $platform ) = @_;
  
  # check to see if it already exists.
  my $run_id = fetch_run_id( $run );
  return $run_id if ( $run_id );

  my $q = "REPLACE INTO run (name, platform) VALUES (?,?)";
  my $sth = $dbi->prepare($q);
  $sth->execute( $run, $platform ) || die "$DBI::errstr";

  return $dbi->last_insert_id(undef, undef, qw( run rid ));
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_runs {
  
  my $q = "SELECT * FROM run";
  my $sth = $dbi->prepare($q);
  $sth->execute( ) || die "$DBI::errstr";

  my @res;
  while (  my @line =  $sth->fetchrow_array()) {
    push @res, [@line];
  }
  
  return @res;
}



# 
# 
# 
# Kim Brugger (24 Jun 2011)
sub fetch_run_name {
  my ( $rid ) = @_;
  my $q = "SELECT name FROM run where rid = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $rid ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return $line[0];
}


# 
# 
# 
# Kim Brugger (24 Jun 2011)
sub fetch_run_id {
  my ( $run ) = @_;
  my $q = "SELECT rid FROM run where name = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $run ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return $line[0];
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

  my $q = "REPLACE INTO project (name) VALUES (?)";
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
# Kim Brugger (23 Jun 2011)
sub fetch_project_name {
  my ( $pid ) = @_;
  my $q = "SELECT name FROM project where pid = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $pid ) || die "$DBI::errstr";
  my @line =  $sth->fetchrow_array();
  return $line[0] || undef;
}



# 
# 
# 
# Kim Brugger (23 Jun 2011)
sub fetch_projects {
  my $q = "SELECT * FROM project";
  my $sth = $dbi->prepare($q);
  $sth->execute( ) || die "$DBI::errstr";
  my @data;
  while (my @line =  $sth->fetchrow_array()) {
    push @data, [@line];
  }

  return @data;
}



# 
# 
# 
# Kim Brugger (20 Jul 2011)
sub fetch_files_from_project {
  my ($pid) = @_;

  my $q = "select s.*, f.fid, f.name from file f, sample s where s.pid= ? and s.sid = f.sid ";
  my $sth = $dbi->prepare($q);
  $sth->execute( $pid ) || die "$DBI::errstr";
  my @data;
  while (my @line =  $sth->fetchrow_array()) {
    push @data, [@line];
  }

  return @data;
}

# 
# 
# 
# Kim Brugger (20 Jul 2011)
sub fetch_files_from_sample {
  my ($sid) = @_;

  my $q = "select f.fid, f.name, r.* from file f, run r where f.sid= ? and r.rid = f.rid";
  my $sth = $dbi->prepare($q);
  $sth->execute( $sid ) || die "$DBI::errstr";
  my @data;
  while (my @line =  $sth->fetchrow_array()) {
    push @data, [@line];
  }

  return @data;
}





1;



__END__




