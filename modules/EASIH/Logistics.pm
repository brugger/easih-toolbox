package EASIH::Logistics;
# 
# Logistics function for checking file names etc. Will most likely be renamed to something better
# at a later point.
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
  $dbname ||= "easih_logistics_dev";
  $dbhost ||= "mgpc17";
#  $dbhost = "localhost";

  $dbi = DBI->connect("DBI:mysql:$dbname:$dbhost", 'easih_admin','easih') || die "Could not connect to database: $DBI::errstr";
}



# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub add_file_offload {
  my ($run_folder, $infile, $outfile) = @_;
  my $sth = $dbi->prepare("INSERT INTO offloads (runfolder, in_file, out_file) VALUES (?,?,?) ");
  $sth->execute( $run_folder, $infile, $outfile );
}


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_files_from_rundir {
  my ($run_folder) = @_;
  my $sth = $dbi->prepare("SELECT out_file FROM offloads WHERE runfolder = ?");
  $sth->execute( $run_folder );
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
sub fetch_offloaded_files {

  my $sth = $dbi->prepare("SELECT runfolder, out_file FROM offloads");
  $sth->execute( ) || die "$DBI::errstr";
  my @res;
  while ( my @line =  $sth->fetchrow_array()) {
    push @res, [@line];
  }

  return @res;
}



# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub runfolder_log {
  my ($run_folder) = @_;

  my $entries = fetch_run_folder_entries($run_folder);
  my @files    = fetch_files_from_rundir($run_folder);

  my $printstring;
  
  $printstring .= "========================================\n";
  $printstring .= "  $run_folder log\n";
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
sub add_run_folder_status {
#  return;
  my ($run_folder, $platform, $status) = @_;
  my $timestamp = Time::HiRes::gettimeofday()*100000;
  my $sth = $dbi->prepare("INSERT INTO runs (runfolder, platform, stamp, status) VALUES (?,?,?,?) ");
  $sth->execute( $run_folder, $platform, $timestamp, $status );
}




# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_run_folder_entries {
  my ($run_folder) = @_;
  
  my $q = "SELECT * FROM runs where runfolder = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $run_folder );
  
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
sub fetch_latest_run_folder_status {
  my ($run_folder) = @_;

  my $entries = fetch_run_folder_entries( $run_folder);
  my $status;
  $status =  $$entries[-1][1] if ($$entries[-1]);
  return $status;
}


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub fetch_platform {
  my ($run_folder) = @_;

  my $q = "SELECT platform from  runs where runfolder = ?";
  my $sth = $dbi->prepare($q);
  $sth->execute( $run_folder );

  my $platform = undef;
  my $results = $sth->fetchrow_arrayref();
  $platform = $$results[0] if ($results);
  return $platform;
}

  
# 
# 
# 
# Kim Brugger (22 Jun 2011)
sub fetch_run_folder_failure {
  my ($run_folder) = @_;

  my $entries = fetch_run_folder_entries( $run_folder);
  my $status;
  $status =  $$entries[-2][1] if ($$entries[-2]);
  return $status;
}



BEGIN{ 
  Connect();
}

1;



