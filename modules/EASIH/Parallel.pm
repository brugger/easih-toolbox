package EASIH::Parallel;
# 
# For running things in parallel...
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

# Gives an error if it is not the master branch + gives access to version information.
use EASIH::Git;

use POSIX ':sys_wait_h';


my @fhs;
my @jobs;
my $MAX_NODES = 8;



# 
# 
# 
# Kim Brugger (20 Sep 2011)
sub max_nodes {
  my ($max_nodes) = @_;

  $MAX_NODES = $max_nodes if (defined  $max_nodes );
  
  return $MAX_NODES;
}


# 
# 
# 
# Kim Brugger (20 Sep 2011)
sub job_push {
  my (@params) = @_;
  
  push @jobs, [@params];

}



# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub run_serial {

  my $output;

  while( my $param = shift @jobs) {
    my $function = shift @$param;
    &$function( @$param);
  }
}


# 
# 
# 
# Kim Brugger (04 Jan 2011)
sub run_parallel {  

  my $output;

  my @cpids     = ();
  
  my $done = 0; # to track the number of files handled

  my $running_nodes = 0;
  my $total = 0;

  while (1) {
  FREE_NODE:
    if ($running_nodes <  $MAX_NODES) {
      my $param = shift @jobs;

      last if ( !$param);

      my $function = shift @$param;
      
      my $cpid = create_child($function, @$param);
      $running_nodes++;
      $total++;
      push @cpids, $cpid;
    }
    else {
      # loop through the nodes to see when one becomes available ...
      while ($done < $total) {
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
	sleep 1;
	
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
  }

  foreach my $fh (@fhs ) {
    next if (! $fh);
    $output .= join("", <$fh>);
    close $fh;
  }
  
  return $output;
}



sub create_child {
  my ($command, @params) = @_;

  my $pid;
  if ($pid = open($fhs[$params[0]], "-|")) {
    ;
  } 
  else {
    die "cannot fork: $!" unless defined $pid;
    &$command(@params);
    exit;
  }
  
  return \$pid;
}







1;



