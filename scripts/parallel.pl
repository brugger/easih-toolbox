#!/usr/bin/perl -w
# 
# Program that can take a list programs to be run and run them 
# parallel on descent computers.
# 
# Nicked and adapted from old cph code.
#
# Kim Brugger (30 Nov 2010), contact: kim.brugger@easih.ac.uk

use strict;
use POSIX ':sys_wait_h';
use POSIX  'tmpnam';
use Getopt::Std;

my %opts;
getopts("c:ns:", \%opts);



my $cpus = nr_of_cpus();
my $MAX_NODES = $opts{c} || $cpus;
my $infile = shift || die "USAGE $0 -c[pus to use, (runs on $cpus cpus by default)] COMMAND-INFILE \n";

my $sleep_time = $opts{s} || 4;

# Store the user specifed values in more readable named variables.
my $INFILE    = $infile;

my @cpids = ();
my (@outfiles, @errfiles);

# Splits the infile up so there is one entry pr. file.

my @running = ("|","/","-","\\", "|", "/", "-");
my $running_counter = 0;

$|  = 1; #Flush all file handles, so we can monitor the progress.
my $done = 0; # to track the number of files handled

my $total = 0;
my $running_nodes = 0;

my @commands;
open (my $in, $infile) || die "Could not open '$infile': $!\n";
while( <$in>) {
  next if (/^\s*\z/);
  next if (/^\#/);
  chomp;
  push @commands, $_;
  $total++;
}  

print "       Running the '$INFILE' infile on '$MAX_NODES' CPUs\n";
print "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n";
print "     Running:      Total :           Done:               Percent:\n";
printf( "         %1s       %5d           %7d                   %3.2f  \r",
	$running[$running_counter++ % 7], $total, $done, $done*100/$total);


while ($_ = shift @commands ) {

 FREE_NODE:
  my $load = load();
  if ( $load < $cpus && $running_nodes <  $MAX_NODES) {
    
      my $command = "$_";
#    $command .= "> /dev/null 2>/dev/null";
      
#    $total++;

      if ( $opts{n}) {
	print "$command\n";
	next;
      }
      
      my $cpid = create_child($command);
      $running_nodes++;
      push @cpids, $cpid;
  }
  else {
    # loop through the nodes to see when one becomes available ...
    while ($running_nodes) {
      for (my $i = 0; $i <@cpids; $i++) {
	next if ($cpids[$i] == -10);
	
	my $cpid = $cpids[$i];
	if (!waitpid($$cpid, WNOHANG)) {
#	  print "Waiting for ($$cpid)\n";
	}
	elsif ($$cpid != -10) {
	  $done++;
	  $cpids[$i] = -10;
	  $running_nodes--;
	}
      }
      sleep $sleep_time;
      printf( "         %1s       %5d           %7d                   %3.2f  \r",
              $running[$running_counter++ % 7], $total, $done, $done*100/$total);
      last if ($running_nodes < $MAX_NODES);
    }
    goto FREE_NODE;
  }

  last if ($running_nodes == 0);
}


while ($done < $total) {
  for (my $i = 0; $i <@cpids; $i++) {
    next if ($cpids[$i] == -10);
    
    my $cpid = $cpids[$i];
    if (!waitpid($$cpid, WNOHANG)) {
      ;
    }
    elsif ($$cpid != -10) {
      $done++;
      $cpids[$i] = -10;
    }
  }
  printf( "         %1s       %5d           %7d                   %3.2f \r",
	  $running[$running_counter++ % 7], $total, $done, $done*100/$total);
  sleep $sleep_time;
}

print "All done.\n";

sub create_child {
  my ($command) = @_;

  my $pid;
  if ($pid = fork) {
    ;
  } 
  else {
    die "cannot fork: $!" unless defined $pid;

    # if the process crashes, run it again, with a limit of 2 times...
    my $limit = 2;
    while (system($command)){
      last if ($limit-- == 0);
    };
    exit;
  }
  
  return \$pid;
}




# 
# 
# 
# Kim Brugger (13 Jan 2011)
sub load {
  
  my $uptime = `uptime`;
  
  $uptime =~ /load average: (\d+.\d+)/;
  return $1 || undef;
}



# 
# 
# 
# Kim Brugger (13 Jan 2011)
sub nr_of_cpus {

  my $cpus = `cat /proc/cpuinfo | egrep ^proc | wc -l`;
  chomp $cpus;
  return $cpus;
}


