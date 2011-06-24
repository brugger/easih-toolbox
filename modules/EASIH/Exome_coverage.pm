package EASIH::Exome_coverage;
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

my $ref_id;

my ($sth_insert_sample, $sth_fetch_sample_id, $sth_fetch_sample_reads);
my ($sth_insert_depth, $sth_fetch_depths);
my $sth_insert_depths;

#
#
#
# Kim Brugger (16 Feb 2011)
sub add_depth {
  my ( $sample_id, $chr, $pos, $depth) = @_;
  die "needs more parameters: $sample_id, $chr, $pos, $depth \n" if ( !$sample_id || !$chr || !$pos || !$depth);

  $sth_insert_depth->execute($sample_id, $chr, $pos, $depth) || die "Insert error\n";;
}



#
#
#
# Kim Brugger (16 Feb 2011)
sub add_depths {
  my ( $sample_id, $chr, $start, $end, $depths) = @_;
  die "needs more parameters: $sample_id, $chr, $start, $end, $depths \n" if ( !$sample_id || !$chr || !$start || !$end || !$depths);

  $depths = join(",", @$depths);

  $sth_insert_depths->execute($sample_id, $chr, $start, $end, $depths) || die "Insert error\n";;
}


# 
# 
# 
# Kim Brugger (25 Feb 2011)
sub fetch_depths {
  my ( $sample_id, $chr, $start, $end) = @_;

  use Data::Dumper;

  my @depths;
  for(my $i=0;$i< $end-$start+1;$i++) {
    push @depths, [$start+$i, 0];
  }

  $sth_fetch_depths->execute( $sample_id,  $chr, $start, $end, $start, $end, $start, $end, $start, $end );

  while (my ($dbstart, $dbend, $dbdepths) =  $sth_fetch_depths->fetchrow_array()) {
    my @dbdepths = split(",", $dbdepths);
    for(my$i=0; $i < $dbend-$dbstart + 1; $i++) {
      next if ($dbstart + $i < $start);
      last if ($dbstart + $i > $end );
      $depths[$dbstart - $start + $i] = [$dbstart+$i, $dbdepths[$i-1]];
    }
  }




  return \@depths;

}



#
#
#
# Kim Brugger (16 Feb 2011)
sub add_sample {
  my ( $name, $reads_mapped, $platform) = @_;
  die "needs more parameters: $name, $reads_mapped, $platform, $ref_id \n" if ( ! $name || !$reads_mapped || !$platform || !$ref_id);

  my $platform_id = set_platform($platform);
  $sth_insert_sample->execute($name, $reads_mapped, $ref_id, $platform_id) || die "Insert error\n";;
  my $sample_id = $sth_insert_sample->{mysql_insertid};

  return $sample_id || undef;
}


# 
# 
# 
# Kim Brugger (16 Feb 2011)
sub fetch_sample_id {
  my ($name) = @_;

  $sth_fetch_sample_id->execute($name);
  my @entries =  $sth_fetch_sample_id->fetchrow_array();
  if (@entries) {
    return $entries[0];
  }
  
  return undef;
}



# 
# 
# 
# Kim Brugger (16 Feb 2011)
sub fetch_sample_reads {
  my ($sample_id) = @_;

  $sth_fetch_sample_reads->execute($sample_id);
  my @entries =  $sth_fetch_sample_reads->fetchrow_array();
  if (@entries) {
    return $entries[0];
  }
  
  return undef;
}




# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub set_platform {
  my ($description) = @_;

  my $sth = $dbi->prepare("select platform_id from platform where description='$description';");
  $sth->execute();
  my @entries =  $sth->fetchrow_array();
  if (@entries) {
    return $entries[0];
  }

  $sth = $dbi->prepare("INSERT INTO platform (description) VALUES ('$description');");
  $sth->execute() || die "Insert failed\n";
  return ($sth->{mysql_insertid});
}


# 
# 
# 
# Kim Brugger (08 Dec 2010)
sub set_ref_id {
  my ($name) = @_;

  my $sth = $dbi->prepare("select ref_id from reference where name='$name';");
  $sth->execute();
  my @entries =  $sth->fetchrow_array();
  if (@entries) {
    $ref_id = $entries[0];
    return $entries[0];
  }

  $sth = $dbi->prepare("INSERT INTO reference (name) VALUES ('$name');");
  $sth->execute() || die "Insert failed\n";
  $ref_id = $sth->{mysql_insertid};
  return ($sth->{mysql_insertid});
}




# 
# 
# 
# Kim Brugger (16 Feb 2011)
sub connect {
  my ($dbase, $user, $password ) = @_;
  
  $dbase ||= "exome_coverage";
  $user ||= 'easih_ro';
  $password ||= '';

  $dbi = DBI->connect("DBI:mysql:$dbase:mgpc17", $user, $password) || die "Could not connect to database: $DBI::errstr";
  $sth_insert_sample      = $dbi->prepare("INSERT INTO sample (name, mapped_reads, ref_id, platform_id) VALUES (?,?,?,?);");
  $sth_fetch_sample_id    = $dbi->prepare("SELECT sample_id FROM sample where name=?;");

  $sth_insert_depth       = $dbi->prepare("INSERT INTO depths (sample_id, chr, pos, depth) VALUES (?,?,?,?);");
  $sth_insert_depths      = $dbi->prepare("INSERT INTO depths (sample_id, chr, start, end, depths) VALUES (?,?,?,?,?);");

  $sth_fetch_depths       = $dbi->prepare("SELECT start, end, depths FROM depths WHERE  sample_id=? AND chr=? AND ((start <= ? AND end >= ?) OR (start >= ? AND end <= ?) OR (start >= ? AND start <= ?)  OR (end >= ? AND end <= ?));");
  $sth_fetch_sample_reads = $dbi->prepare("SELECT mapped_reads FROM sample WHERE sample_id=?;");
}

1;
