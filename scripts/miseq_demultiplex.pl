#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (18 Jul 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my %files;
my %out_fhs;
input_files();
my $sample_sheet = sample_sheet();


my ($r1_fh, $r2_fh, $i1_fh, $i2_fh);

# print shift @{$files{ R1 }}, "\n";
# print shift @{$files{ R2 }}, "\n";
# print shift @{$files{ I1 }}, "\n";
# print shift @{$files{ I2 }}, "\n";


open_file(\$r1_fh, shift @{$files{ R1 }});
open_file(\$r2_fh, shift @{$files{ R2 }});
open_file(\$i1_fh, shift @{$files{ I1 }});
open_file(\$i2_fh, shift @{$files{ I2 }});

while (1) {
  
  my ($headr1,$seqr1,$strr1,$qualr1) = next_entry( $r1_fh );
  if ( !$headr1 && @{$files{ R1 }}) {
    close($r1_fh);
    close($r2_fh);
    close($i1_fh);
    close($i2_fh);
    open_file(\$r1_fh, shift @{$files{ R1 }});
    open_file(\$r2_fh, shift @{$files{ R2 }});
    open_file(\$i1_fh, shift @{$files{ I1 }});
    open_file(\$i2_fh, shift @{$files{ I2 }});

    next;
  }
  elsif ( !$headr1 ) {
    exit;
  }
  my ($headr2,$seqr2,$strr2,$qualr2) = next_entry( $r2_fh );
  my ($headi1,$seqi1,$stri1,$quali1) = next_entry( $i1_fh );
  my ($headi2,$seqi2,$stri2,$quali2) = next_entry( $i2_fh );
  
  if (0) {
    print "$headr1$seqr1$strr1$qualr1";
    print "$headr2$seqr2$strr2$qualr2";

    print "$headi1$seqi1$stri1$quali1";
    print "$headi2$seqi2$stri2$quali2";
  }

  chomp $seqi1;
  chomp $seqi2;

  if ( $$sample_sheet{ $seqi1 }{ $seqi2 }) {
    my ($fh1, $fh2) = fetch_fh($$sample_sheet{ $seqi1 }{ $seqi2 });

    print $fh1 "$headr1$seqr1$strr1$qualr1";
    print $fh2 "$headr2$seqr2$strr2$qualr2";
  }


#  exit;
}



# 
# 
# 
# Kim Brugger (18 Jul 2012)
sub fetch_fh {
  my ($name ) = @_;

  return ($out_fhs{ $name }{1}, $out_fhs{ $name }{2}) 
      if ( $out_fhs{ $name } && $out_fhs{ $name }{1});
  open ( $out_fhs{ $name }{1}, "> $name.1.fq") || die "Could not write to $name.1.fq\n";
  open ( $out_fhs{ $name }{2}, "> $name.2.fq") || die "Could not write to $name.1.fq\n";

  return ($out_fhs{ $name }{1}, $out_fhs{ $name }{2});
  
}


# 
# 
# 
# Kim Brugger (18 Jul 2012)
sub next_entry {
  my( $fh ) = @_;

  my $head = <$fh>;
  
  return (undef) if ( ! $head );


  my $seq  = <$fh>;
  my $str  = <$fh>;
  my $qual = <$fh>;
  
  return ($head,$seq,$str,$qual);
}



# 
# 
# 
# Kim Brugger (18 Jul 2012)
sub open_file {
  my ($fh, $file) = @_;
  open( $$fh, "gunzip -c $file | ") || die "Could not open file '$file': $!\n";
  
  
}




# 
# 
# 
# Kim Brugger (18 Jul 2012)
sub sample_sheet {
  my $in = shift;
  $in ||= "SampleSheet.csv";


  my %res;

  open (my $ss, "$in") || die "Could not open file '$in': $!\n";

  while (<$ss>) {
    chomp;
    s/\r//g;
    s/\n//g;

    my ($name, $bc1, $bc2) = split(",");
    
    $res{$bc1}{$bc2} = $name;

  }  
  return \%res;
}


# 
# 
# 
# Kim Brugger (18 Jul 2012)
sub input_files {
  my ($indir) = @_;
  $indir .= "/" if ($indir);
  $indir ||= "";

  foreach my $file ( glob("$indir*gz")) {
    if ( $file =~ /R1/ ) {
      push @{$files{ R1 }}, $file;
    }
    elsif ( $file =~ /R2/ ) {
      push @{$files{ R2 }}, $file;
    }
    elsif ( $file =~ /I1/ ) {
      push @{$files{ I1 }}, $file;
    }
    elsif ( $file =~ /I2/ ) {
      push @{$files{ I2 }}, $file;
    }
    
  }
  
  
  foreach my $key ( keys %files ) {
    @{$files{$key}} = sort @{$files{$key}};
  }
  
}

