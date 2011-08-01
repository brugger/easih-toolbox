package EASIH::Sample;
# 
# Generic function for handling sample names and extract project ID etc from them
# 
# 
# Kim Brugger (30 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub validate_name {
  my ($sample_name) = @_;

  # Takes the old names into consideration as well.
  return 1 if ($sample_name =~ /^[A-Z]\d{6,7}\z/);

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



# 
# 
# 
# Kim Brugger (27 Jun 2011)
sub filename2sampleNproject {
  my ($filename) = @_;

  # remove the path
  $filename =~ s/.*\///;

  # extract the sample name
  my ($sample) = $filename =~ /^([A-Z]\d{6,7})/;
  if ( ! $sample ) {
    print STDERR "$sample does not comply with the sample naming scheme\n";
    return (undef,undef);
  }

  # should _never_ happen as is should have been caught above.
  my ($project) = $sample =~ /^([A-Z]\d{2})/;
  if ( ! $project ) {
    print STDERR "$project does not comply with the sample naming scheme\n";
    return (undef, undef);
  }

  return ($sample, $project);
}

1;


