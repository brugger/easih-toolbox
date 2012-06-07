package EASIH::ChiSquare;

# ChiSquare.pm
#
# merged and extended from two existing perl modules (Statistics::Distributions and Statistics::ChiSquare)
#
#
# Kim Brugger (20 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use vars qw($VERSION @ISA @EXPORT);

require Exporter;
use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; # number of significant digits to be returned

@ISA = qw(Exporter);
@EXPORT = qw(chisquare chisquare_unbalanced);

$VERSION = '0.5';

# assume the expected probability distribution is uniform
sub chisquare {
  my @data = @_;
  @data = @{$data[0]} if @data == 1 and ref($data[0]);
  return "There's no data!" unless @data;
  
  my $degrees_of_freedom = scalar(@data) - 1;
  my ($chisquare, $num_samples, $expected) = (0, 0, 0);
  foreach (@data) { 
    $num_samples += $_ 
  }
  $expected = $num_samples / scalar(@data);
  # return "There's no data!" unless $expected;
  foreach (@data) {
    $chisquare += (($_ - $expected) ** 2) / $expected;
  }
  $chisquare = sprintf("%.5f", $chisquare);
  

  my $chip = chisqrprob ($degrees_of_freedom,$chisquare);
  return "[chi-square: $chisquare, df: $degrees_of_freedom] ==>  P = $chip";
}

sub chisquare_raw {
  my @data = @_;
  @data = @{$data[0]} if @data == 1 and ref($data[0]);
  return "There's no data!" unless @data;
  
  my $degrees_of_freedom = scalar(@data) - 1;
  my ($chisquare, $num_samples, $expected) = (0, 0, 0);
  foreach (@data) { 
    $num_samples += $_ 
  }
  $expected = $num_samples / scalar(@data);
  # return "There's no data!" unless $expected;
  foreach (@data) {
    $chisquare += (($_ - $expected) ** 2) / $expected;
  }
  $chisquare = sprintf("%.5f", $chisquare);
  

  my $chip = chisqrprob ($degrees_of_freedom,$chisquare);
  return ($chisquare, $degrees_of_freedom, $chip);
}


# assume the expected probability distribution is uniform
sub chisquare_unbalanced {
  my ($data, $freq) = @_;
  my @data = @$data;
  my @freq = @$freq;

  @data = @{$data[0]} if @data == 1 and ref($data[0]);
  return "There's no data!" unless @data;
  
  my $degrees_of_freedom = scalar(@data) - 1;
  my ($chisquare, $num_samples) = (0, 0 );
  foreach (@data) { 
    $num_samples += $_ 
  }
  # return "There's no data!" unless $expected;
  for(my $i = 0; $i< @data; $i++) {
    my $count = $data[$i];
    my $expected = $num_samples * $freq[$i];
    $chisquare += (($count - $expected) ** 2) / $expected;
  }
  $chisquare = sprintf("%.5f", $chisquare);
  

  my $chip = chisqrprob ($degrees_of_freedom,$chisquare);
  return "[chi-square: $chisquare, df: $degrees_of_freedom] ==>  P = $chip";
}



   
sub chisqrdistr { # Percentage points  X^2(x^2,n)
	my ($n, $p) = @_;
	if ($n <= 0 || abs($n) - abs(int($n)) != 0) {
		die "Invalid n: $n\n"; # degree of freedom
	}
	if ($p <= 0 || $p > 1) {
		die "Invalid p: $p\n"; 
	}
	return precision_string(_subchisqr($n, $p));
}

sub chisqrprob { # Upper probability   X^2(x^2,n)
	my ($n,$x) = @_;
	if (($n <= 0) || ((abs($n) - (abs(int($n)))) != 0)) {
		die "Invalid n: $n\n"; # degree of freedom
	}
	return precision_string(_subchisqrprob($n, $x));
}



sub _subchisqrprob {
	my ($n,$x) = @_;
	my $p;

	if ($x <= 0) {
		$p = 1;
	} elsif ($n > 100) {
		$p = _subuprob((($x / $n) ** (1/3)
				- (1 - 2/9/$n)) / sqrt(2/9/$n));
	} elsif ($x > 400) {
		$p = 0;
	} else {   
		my ($a, $i, $i1);
		if (($n % 2) != 0) {
			$p = 2 * _subuprob(sqrt($x));
			$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
			$i1 = 1;
		} else {
			$p = $a = exp(-$x/2);
			$i1 = 2;
		}

		for ($i = $i1; $i <= ($n-2); $i += 2) {
			$a *= $x / $i;
			$p += $a;
		}
	}
	return $p;
}


sub _subuprob {
	my ($x) = @_;
	my $p = 0; # if ($absx > 100)
	my $absx = abs($x);

	if ($absx < 1.9) {
		$p = (1 +
			$absx * (.049867347
			  + $absx * (.0211410061
			  	+ $absx * (.0032776263
				  + $absx * (.0000380036
					+ $absx * (.0000488906
					  + $absx * .000005383)))))) ** -16/2;
	} elsif ($absx <= 100) {
		for (my $i = 18; $i >= 1; $i--) {
			$p = $i / ($absx + $p);
		}
		$p = exp(-.5 * $absx * $absx) 
			/ sqrt(2 * PI) / ($absx + $p);
	}

	$p = 1 - $p if ($x<0);
	return $p;
}


sub _subchisqr {
	my ($n, $p) = @_;
	my $x;

	if (($p > 1) || ($p <= 0)) {
		die "Invalid p: $p\n";
	} elsif ($p == 1){
		$x = 0;
	} elsif ($n == 1) {
		$x = _subu($p / 2) ** 2;
	} elsif ($n == 2) {
		$x = -2 * log($p);
	} else {
		my $u = _subu($p);
		my $u2 = $u * $u;

		$x = max(0, $n + sqrt(2 * $n) * $u 
			+ 2/3 * ($u2 - 1)
			+ $u * ($u2 - 7) / 9 / sqrt(2 * $n)
			- 2/405 / $n * ($u2 * (3 *$u2 + 7) - 16));

		if ($n <= 100) {
			my ($x0, $p1, $z);
			do {
				$x0 = $x;
				if ($x < 0) {
					$p1 = 1;
				} elsif ($n>100) {
					$p1 = _subuprob((($x / $n)**(1/3) - (1 - 2/9/$n))
						/ sqrt(2/9/$n));
				} elsif ($x>400) {
					$p1 = 0;
				} else {
					my ($i0, $a);
					if (($n % 2) != 0) {
						$p1 = 2 * _subuprob(sqrt($x));
						$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
						$i0 = 1;
					} else {
						$p1 = $a = exp(-$x/2);
						$i0 = 2;
					}

					for (my $i = $i0; $i <= $n-2; $i += 2) {
						$a *= $x / $i;
						$p1 += $a;
					}
				}
				$z = exp((($n-1) * log($x/$n) - log(4*PI*$x) 
					+ $n - $x - 1/$n/6) / 2);
				$x += ($p1 - $p) / $z;
				$x = sprintf("%.5f", $x);
			} while (($n < 31) && (abs($x0 - $x) > 1e-4));
		}
	}
	return $x;
}

sub log10 {
	my $n = shift;
	return log($n) / log(10);
}
 
sub max {
	my $max = shift;
	my $next;
	while (@_) {
		$next = shift;
		$max = $next if ($next > $max);
	}	
	return $max;
}

sub min {
	my $min = shift;
	my $next;
	while (@_) {
		$next = shift;
		$min = $next if ($next < $min);
	}	
	return $min;
}

sub precision {
	my ($x) = @_;
	return abs int(log10(abs $x) - SIGNIFICANT);
}

sub precision_string {
	my ($x) = @_;
	if ($x) {
		return sprintf "%." . precision($x) . "f", $x;
	} else {
		return "0";
	}
}


1;
