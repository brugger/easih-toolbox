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

  print "========================================\n";
  print "  $run_folder log\n";
  print "========================================\n";

  use POSIX qw( strftime );
  foreach my $entry ( @$entries ) {
    my $microsec = $$entry[0] % 100000;
    $$entry[0] /= 100000;
    $$entry[0] = strftime("%Y/%m/%d %H:%M:%S:$microsec", localtime( $$entry[0] ));
    print "$$entry[0]\t$$entry[1]\n";
  }

  print "========================================\n";
  foreach my $file ( @files ) {
    print "$file\n";
  }
  print "========================================\n";



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
sub validate_sample_name {
  my ($sample_name) = @_;

  # Takes the old names into consideration as well.
  return 1 if ($sample_name =~ /[A-Z]\d{6,7}/);

  return 0;
}



# 
# Finds the filename a sample should go to. Will check for std naming + if previous files exsists.
# outdir overrides the default /data/<PROJECTID>/raw folder.
#
# It will furthermore check for permissions and create missing directories.
# 
# Kim Brugger (13 Jun 2011)
sub sample2outfilename {
  my ($sample_name, $outdir) = @_;

  chomp($sample_name);

  # replace the torrent fastq with a fq to perserve consistent naming
  $sample_name =~ s/fastq/fq/;
  
  my ($sample, $postfix) = $sample_name =~ /^([A-Z]\d{6,7})(.*)/;
  
  if ( $outdir ) {
    if ( -e "$outdir" && ! -d "$outdir") {
      return (undef,  "$outdir is not a directory\n");
#      die "$outdir is not a directory\n";
    }
    elsif ( -e "$outdir" && ! -w "$outdir") {
      return (undef,  "$outdir is not writeable\n");
#      die "$outdir is not writeable\n";
    }
    if ( $outdir && ! -d $outdir ) {
      if (!system "mkdir -p $outdir") {
	return (undef,  "Could not create directory '$outdir': $!\n");
      }
    }
  }

  my $root_dir = "/data/";
  my $project = substr($sample, 0, 3);
  my $file = "$root_dir/$project/raw/$sample";

  if ( $outdir ) {
    $root_dir = $outdir;
    $file     = "$outdir/$sample";
  }
  else {
    if ( -e "$root_dir/$project" && ! -d "$root_dir/$project") {
      return (undef, "$root_dir/$project is not a directory\n");
    }
    elsif ( -e "$root_dir/$project" && ! -w "$root_dir/$project") {
      return (undef, "$root_dir/$project is not writeable\n");
    }
    
    $root_dir .= "$project/raw/";
    system "mkdir -p $root_dir" if ( ! -d "$root_dir" );
  }
  
  
  my @files = `find $root_dir | grep $sample `;
  my $version = 0;
  if ( @files ) {
    while ( $_ = pop @files ) {
      chomp;
      if (  /[A-Z]\d{6,7}\_(\d+)\.\d+$postfix/ || /[A-Z]\d{6,7}\_(\d+)$postfix/) {
	$version = $1+1 if ($version < $1+1 );
      }
      elsif ( $version == 0 && ( /[A-Z]\d{6,7}\.\d+$postfix/  || /[A-Z]\d{6,7}$postfix/)) {
	$version = 1;
      }
    }
    
    $file = "$root_dir/$sample\_$version";
  }

  return ("$file$postfix", undef);
}



BEGIN{ 
  Connect();
}

1;



