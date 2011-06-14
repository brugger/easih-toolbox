#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (12 Jan 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my @pre;

my (@fields);
my %field_names;

my %effects  = ('ESSENTIAL_SPLICE_SITE'  => 10, 
		'STOP_GAINED'            => 10, 
		'STOP_LOST'              => 10, 
		'COMPLEX_INDEL'          => 10,
		'FRAMESHIFT_CODING'      => 10, 
		'NON_SYNONYMOUS_CODING'  => 10, 
		'SPLICE_SITE'            =>  8, 
		'PARTIAL_CODON'          =>  8,
		'SYNONYMOUS_CODING'      =>  1, 
		'REGULATORY_REGION'      =>  3, 
		'WITHIN_MATURE_miRNA'    =>  3, 
		'5PRIME_UTR'             =>  1,
		'3PRIME_UTR'             =>  1, 
		'UTR'                    =>  1, 
		'INTRONIC'               =>  1, 
		'NMD_TRANSCRIPT'         =>  1, 
		'WITHIN_NON_CODING_GENE' =>  5, 
		'UPSTREAM'               =>  2,
		'DOWNSTREAM'             =>  2, 
		'HGMD_MUTATION'          => 10, 
		'INTERGENIC'             =>  0,
		'NO_CONSEQUENCE'         =>  0, 
    );


while(<>) {
  chomp;
  my @f = split("\t");

  next if (/#/);
  if( /^position/) {

    my @line = ("chr", "pos");
    
    for(my $i=0; $i< @f; $i++ ) {	
      $field_names{ $f[$i] } = $i;
      
      next if ( ! $f[$i]);

      if ($f[$i] eq "change" ||  
	  $f[$i] eq "filter" ||
	  $f[$i] eq "score" ||   
	  $f[$i] eq "depth" ||  
	  $f[$i] eq "type" ||                                            
	  $f[$i] eq "gene" ||
	  $f[$i] eq "region" ||  
	  $f[$i] eq "codon pos" ||       
	  $f[$i] eq "AA change" ||       
	  $f[$i] eq "Grantham score" || 
	  $f[$i] eq "external ref" ||    
	  $f[$i] eq "dbsnp flags" ||
	  $f[$i] eq "pfam") {
	push @line, $f[$i];
	push @fields, $i;
      }
    }
    
    print join("\t", @line) ."\n";
    next;
  }


  if ( @pre && $f[0] ne $pre[0] ) {
    print_line(@pre);
    @pre = ();
  }
  elsif ( @pre ) {

    my $old_region = $pre[$field_names{'region'}];
    my $new_region = $f[$field_names{'region'}];

    @pre = @f if ( $effects{ $new_region} > $effects{ $old_region });
  }
  else {
    @pre = @f;
  }


}

print_line(@pre);




# 
# 
# 
# Kim Brugger (12 Jan 2011)
sub print_line {
  my (@v) = @_;

  my $region = $v[$field_names{'region'}];
  return if (! $region );
  return if ($effects{$region} < 5);

  my ($chr, $pos) = split(":", $v[0]);

  $v[ $field_names{gene} ] =~ s/^(.*?)\/.*/$1/;

  my @line = ($chr, $pos);

  map { push @line, $v[$_] || "" } @fields;
  print join("\t", @line) ."\n";  

}

