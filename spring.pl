#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use File::Temp;

use Getopt::Std;

my %job_stats;

my %opts;
getopts('w:i:o:p:n:', \%opts);

my $infile  = $opts{'i'} || "1_1_1.fastq" || usage();
my $outfile = $opts{'o'} || "1_1_1.fastq.bam" || usage();
my $prefix  = $opts{'p'} || usage();
my $split   = $opts{'n'} || 10000 || 30000000;

my $fq_split = '/home/kb468/bin/fastq_split.pl';
my $bwa      = '/home/kb468/bin/bwa';
my $samtools = '/home/kb468/bin/samtools';
my $cwd      = $opts{'w'} || `pwd`;
chomp($cwd);

my @analysis = ({logic_name => 'fastq-split',
		 command    => "$fq_split -e $split $infile ",
		 hpc_param  => "-NEP-fqs -l mem=500mb,walltime=00:05:00"},
		{logic_name => 'BWA-mapping',
		 command    => "$bwa aln ",
		 hpc_param  => "-NEP-fqs -l mem=1000mb,walltime=00:10:00"},
		{logic_name => 'BWA-samse',
		 command    => "$bwa samse $prefix ",
		 hpc_param  => "-NEP-fqs -l mem=1000mb,walltime=00:10:00"},
		{logic_name => 'SAM2BAM',
		 command    => "$samtools view -b -S ",
		 hpc_param  => "-NEP-fqs -l mem=500mb,walltime=00:10:00"},
		{logic_name => 'BWA-merge',
		 command    => "$samtools  ",
		 hpc_param  => "-NEP-fqs -l mem=500mb,walltime=00:10:00"},
		{logic_name => 'samtools-sort',
		 command    => "$samtools sort -m 2048000000",
		 hpc_param  => "-NEP-fqs -l mem=1000mb,walltime=00:10:00"},
		{logic_name => 'bam-rename',
		 command    => "mv "},
		{logic_name => 'samtools-index',
		 command    => "$samtools index $outfile",
		 hpc_param  => "-NEP-fqs -l mem=1000mb,walltime=00:10:00"},
		);
		     

my @jobs;
my @delete_files;

my @next_inputs;

while (my $next_task = shift  @analysis ) {
  print "Running $$next_task{ logic_name }' \n";

  if ( $$next_task{ logic_name } eq 'fastq-split') {
    my $tmp_file = tmp_file();
    my $cmd = $$next_task{command } . " > $tmp_file";

    my $job_id = submit_job($cmd, $$next_task{hpc_param});
    
    wait_jobs( $job_id );
    @next_inputs = ();
    open (my $tfile, $tmp_file) || die "Could not open '$tmp_file':$1\n";
    while(<$tfile>) {
      chomp;
      push @next_inputs, $_;
    }
    close ($tmp_file);
    print "Next job_ids = @next_inputs\n";
  }
  elsif ($$next_task{ logic_name } eq 'BWA-mapping') {
    my @inputs = @next_inputs;
    @next_inputs = ();

    my @tmp_job_ids;
    foreach my $input ( @inputs ) {
      my $tmp_file = tmp_file();
      my $cmd = $$next_task{command } . " -f $tmp_file $prefix $input ";
      push @next_inputs, "$tmp_file $input";

      my $job_id = submit_job($cmd, $$next_task{hpc_param});
      push @tmp_job_ids, $job_id;
    }

    wait_jobs( @tmp_job_ids );

    print "Next job_ids = @next_inputs\n";
  }
  elsif ($$next_task{ logic_name } eq 'BWA-samse')  {
    my @inputs = @next_inputs;
    @next_inputs = ();
    my @tmp_job_ids;
    foreach my $input ( @inputs ) {
      my $tmp_file = tmp_file(".sam");
      my $cmd = $$next_task{command } . " $input > $tmp_file ";
      push @next_inputs, $tmp_file;

      my $job_id = submit_job($cmd, $$next_task{hpc_param});
      push @tmp_job_ids, $job_id;
    }

    wait_jobs( @tmp_job_ids );

    print "Next job_ids = @next_inputs\n";
  }
  elsif ($$next_task{ logic_name } eq 'SAM2BAM' ) {
    my @inputs = @next_inputs;
    @next_inputs = ();
    my @tmp_job_ids;
    foreach my $input ( @inputs ) {
      my $tmp_file = tmp_file(".bam");
      my $cmd = $$next_task{command } . " $input > $tmp_file ";

      my $job_id = submit_job($cmd, $$next_task{hpc_param});
      push @next_inputs, $tmp_file;
    }
    wait_jobs( @tmp_job_ids );

    print "Next job_ids = @next_inputs\n";
  }
  elsif ($$next_task{ logic_name } eq 'BWA-merge') {
    my @inputs = @next_inputs;
    @next_inputs = ();
    my $tmp_file = tmp_file(".merged.bam");
    my $cmd = $$next_task{command } . "$tmp_file @inputs ";
    push @next_inputs, $tmp_file;
    
    my $job_id = submit_job($cmd, $$next_task{hpc_param});
    wait_jobs( $job_id );

    print "Next job_ids = @next_inputs\n";
  }
  elsif ($$next_task{ logic_name } eq 'samtools-sort' ) {
    my @inputs = @next_inputs;
    @next_inputs = ();
    my $tmp_file = tmp_file(".merged.sorted");
    my $cmd = $$next_task{command } . " @inputs > $tmp_file ";
    push @next_inputs, "$tmp_file.bam";
    
    my $job_id = submit_job($cmd, $$next_task{hpc_param});
    wait_jobs( $job_id );

  }
  elsif ($$next_task{ logic_name } eq 'bam-rename' ) {
    my @inputs = @next_inputs;
    @next_inputs = ();
    my $cmd = $$next_task{command } . " @inputs $outfile ";
  }
  elsif ($$next_task{ logic_name } eq 'samtools-index' ) {
    my $cmd = $$next_task{command };
    my $job_id = submit_job($cmd, $$next_task{hpc_param});
    wait_jobs( $job_id );
  }


}



# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub submit_job {
  my ($cmd, $hpc_params) = @_;
  my ($tmp_fh, $tmp_file) = File::Temp::tempfile(DIR => "./tmp" );

  open (my $qpipe, " | qsub $hpc_params -o q-logs > $tmp_file 2> /dev/null ") || die "Could not open qsub-pipe: $!\n";
  print $qpipe "cd $cwd; $cmd";
  close( $qpipe );

  open (my $tfile, $tmp_file) || die "Could not open '$tmp_file':$1\n";
  my $job_id;
  while(<$tfile>) {
    chomp;
    $job_id = $_;
  }
  close ($tfile);
  system "rm $tmp_file";

  $job_id =~ s/(\d+?)\..*/$1/;
  return $job_id;
}


# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub wait_jobs {
  my (@job_ids) = @_;
  
  my %s2status = ( C =>  "Completed",
                   E =>  "Exiting",
		   F =>  "Failed",
                   H =>  "Halted",
                   Q =>  "Queued",
                   R =>  "Running",
                   T =>  "Moving",
                   W =>  "Waiting",
                   S =>  "Suspend" );


  my %job_hash;
  map { $job_hash{ $_ }{ full_status } = "UNKNOWN"  } @job_ids;


  while (1) {
    
    my ( $done, $running, $waiting, $queued, $other, ) = (0,0,0,0,0);
    foreach my $job ( keys %job_hash ) {
      
      my ($status, $status_hash) = job_stats( $job );

#      print "STATUS $status\n";

      $done++    if ($status eq 'C');
      $running++ if ($status eq 'R');
      $queued++  if ($status eq 'Q');
      $running++ if ($status eq 'W');
      $other++   if ($status ne 'R' && $status ne 'W' && $status ne 'Q' && $status ne 'C');

      if ( $status eq "F") {
	die "$job failed, so exiting the pipeline\n";
      }

      $job_hash{ $job }{ status }      = $status;
      $job_hash{ $job }{ full_status } = $s2status{ $status };
    }

    print "Job tracking stats: D: $done, R: $running, Q: $queued, W: $waiting, O: $other\n";
    last if ( $done == @job_ids);

    sleep(10);

  }      

  return;
}




# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub job_stats {
  my ($job_id)  = @_;


  my %res;
  open (my $qspipe, "qstat -f $job_id 2> /dev/null | ") || die "Could not open 'qstat-pipeline': $!\n";
  my ( $id, $value);
  while(<$qspipe>) {
    chomp;
#    print "$_ \n";
    if (/ *(\w+) = (.*)/ ) {
      $res{$id} = $value if ( $id && $value);
      $id    = $1;
      $value = $2;
    }
    elsif (/\t(.*)/) {
      $value .= $1;
    }
  }
    
#  print Dumper( \%res );

  $res{job_state} = "F" if ( $res{exit_status} && $res{exit_status} != 0);

  return $res{job_state}, \%res;
  
}




# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub tmp_file {
  my ($postfix) = @_;
  $postfix ||= "";
  my ($tmp_fh, $tmp_file) = File::Temp::tempfile(DIR => "./tmp" );

  push @delete_files, "$tmp_file$postfix";

  return "$tmp_file$postfix";
}



# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub delete_tmp_file {

  system "rm @delete_files";
}




# 
# 
# 
# Kim Brugger (22 Apr 2010)
sub usage {
  
  print "not the right usage\n";
#  exit;

}
