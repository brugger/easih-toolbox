#!/usr/bin/perl -w

use strict;
use English qw( -no_match_vars );
use warnings FATAL => 'all';

use Getopt::Long 'GetOptions';
use Data::Dumper;
use Carp;

my $proper_pm_mask      = 2; #0x0002
my $query_unmap_mask    = 4; #0x0004
my $mate_unmap_mask     = 8; #0x0008
my $strand_mask         = 16; #0x0010 (set when reverse)
my $first_in_pair_mask  = 64;
my $second_in_pair_mask = 128;

### these may need to be parameters.
my $max_n_frac          = 0.10;
my $min_seq_hardlimit   = 22;
my $one_edit_per_bps    = 20;
###

my %file_handles;
my %cent_coords;

#local $| = 1;

my $debug = 0;

{
  my $options = option_builder();
  rg_parse($options);
#  print Dumper($options);
  if($options->{'c'}) {
  	centromere_parser($options->{'c'});
  }
  run_parser($options);
}

sub centromere_parser {
	my ($cent_file) = @_;
	open my $CENTS, '<', $cent_file or croak "Could not open file $cent_file";
	while (my $line=<$CENTS>) {
		if($line =~ m/^#/) {
			next;
		}
		my @gff_bits = split /\t/, $line;
		my $cent_chr = $gff_bits[0];
		my $cent_start = $gff_bits[3];
		my $cent_end = $gff_bits[4];
		$cent_coords{$cent_chr} = [$cent_start, $cent_end];
	}
	close $CENTS;
	return;
}

sub rg_parse {
  my ($options) = @_;
  # get the header only
  my $command = 'samtools view -H -F 12 '.$options->{'i'};
  my $pid = open my $PROC, '-|', $command or croak "Could not fork: $OS_ERROR";
  
  while( my $tmp = <$PROC> ) {
    if($tmp !~ /^\@RG/xs) {
      next;
    }
    if(!$options->{'rgs'}) {
      $options->{'rgs'} = {};
    }
    chomp $tmp;
    my @rg_data = split /\t/xms, $tmp;
    # have to find tag each time as merged data may have different tag positions
    my $id_pos = tag_finder('ID', \@rg_data, 1); # no point parsing 0, as just @RG tag
    my $pi_pos = tag_finder('PI', \@rg_data, 1); # no point parsing 0, as just @RG tag
    my (undef, $id) = split /:/xms, $rg_data[$id_pos];
    my (undef, $pi) = split /:/xms, $rg_data[$pi_pos];
    $options->{'rgs'}->{$id} = $pi;
  }
  close $PROC or croak 'Failed to close process: '.$command."\nERROR CODE: ".$CHILD_ERROR;
}

sub run_parser {
  my ($options) = @_;
  #my $command = 'samtools view -F 12 '.$options->{'i'};
  ## this will prevent any pairs with both unmapped
  ## best method is to exclude these from the input by doing
  ## samtools view -F 12 -q 1 -b -o NEW_FILE.bam ORIGINAL.bam
  
  my $command = 'samcat -n ';
  if($options->{'om'}) { # use oldmethod for parsing of bam file
  	$command = 'samtools view ';
  }
  $command .= $options->{'i'};
  
  my ($file_chr, %sw_reads);
  my $pid = open my $PROC, '-|', $command or croak "Could not fork: $OS_ERROR";
  my (@curr_loc, $chr, $pos, $xt_pos);
  my $discarded = 0;
  
  my @pair;
  MAIN: while( my $tmp = <$PROC> ) {
  	if($tmp =~ m/^\@/) {
  		next;
  	}
  	chomp $tmp;
    my @comps = split(/\t/, $tmp);
    @pair = (\@comps);
    
    # this forces a paired set of reads into @pair
    my $paired = 0;
    while($paired == 0) {
    	$tmp = <$PROC>;
    	if(!$tmp) { last MAIN; }
    	chomp $tmp;
    	my @comps2 = split(/\t/, $tmp);
    	if($pair[0]->[0] ne $comps2[0]) {
    		$pair[0] = \@comps2;
    	}
    	else {
    		push @pair, \@comps2;
    		$paired = 1;
    	}
    }
    my ($rec, $rec_2);
    eval {
    	my $proc_res = process_data(\@pair, $options);
    	if($proc_res) {
    		($rec, $rec_2) = @{$proc_res};
    	}
    };
    if($EVAL_ERROR) {
    	my $error =  "\n".$EVAL_ERROR."\n\n";
    	$error .= join "\t", @{$pair[0]};
    	$error .= join "\t", @{$pair[1]};
    	croak $error;
    }
		
		my $CHR_FH;
		if($rec && @{$rec} > 0) {
			$CHR_FH = get_chr_fh($options, $rec->[1], $rec->[2]);
			print $CHR_FH $rec->[0];
			if($rec_2 && @{$rec_2} > 0) {
				$CHR_FH = get_chr_fh($options, $rec_2->[1], $rec_2->[2]);
				print $CHR_FH $rec_2->[0];
			}
			#if($debug) {
				#<STDIN>;
			#}
		}
  }
  # close samtools view
  close $PROC or croak 'Failed to close process: '.$command."\nERROR CODE: ".$CHILD_ERROR;
  
  # close output files
  foreach my $file_keys(sort keys %file_handles) {
  	close $file_handles{$file_keys}
  }
}

sub get_chr_fh {
	my ($options, $chr, $pos) = @_;
	my $file_chr = $chr;
	if($options->{'c'} && $chr ne 'MT') {
		if($pos > $cent_coords{$chr}->[0]) {
			$file_chr .= '-q';
		}
		else {
			$file_chr .= '-p';
		}
	}
	my $CHR_FH;
	if(!$file_handles{$file_chr}) {
		my $file_path = $options->{'o'}.'_'.$file_chr.'.txt';
		open $CHR_FH, '>', $file_path  or croak 'Could not create file: '.$file_path;
		$file_handles{$file_chr} = $CHR_FH;
	}
	else {
		$CHR_FH = $file_handles{$file_chr};
	}
	return $CHR_FH;
}

sub process_data {
  my ($read_pair, $options) = @_;
#  \@curr_loc, $OUTFILE, $chr, $file_chr, $options, $discarded, \%sw_reads
  
  # setup for basic short circuits to skip readpairs
  
  my $chr1 = $read_pair->[0]->[2];
	my $chr2 = $read_pair->[1]->[2];
	
	#if($chr1 ne '16' && $chr2 ne '16') {
	#	return;
	#}
	
	my $comps1 = $read_pair->[0];
	my $comps2 = $read_pair->[1];
	
	my $read_name = $comps1->[0];
	
	
  my $cigar1 = $comps1->[5];
  my $cigar2 = $comps2->[5];
  
  debug(join "\t", @{$comps1});
  debug(join "\t", @{$comps2});
  
  if($cigar1 eq '*' && $cigar2 eq '*') { # both reads unmapped
  	debug($read_name.' Both reads unmapped');
    return;
  }
  
  # anchor read needs to be uniquely mapped
  my ($xt_1, $xt_2) = (0,0);
  my ($best1, $best2) = (0,0);
  my ($subopt1, $subopt2) = (0,0);
  if($cigar1 ne '*') {
		(undef, undef, $xt_1) = split /:/xms, $comps1->[tag_finder('XT', $comps1, 11)];
		if($xt_1 ne 'M' && $xt_1 ne 'N') {
    	(undef, undef, $best1) = split /:/, $comps1->[tag_finder('X0', $comps1, 11)];
    	my $x1_pos = tag_finder('X1', $comps1, 11);
    	if($x1_pos) {
	    	(undef, undef, $subopt1) = split /:/, $comps1->[$x1_pos];
    	}
    }
  }
  if($cigar2 ne '*') {
    (undef, undef, $xt_2) = split /:/xms, $comps2->[tag_finder('XT', $comps2, 11)];
    if($xt_2 ne 'M' && $xt_2 ne 'N') {
    	(undef, undef, $best2) = split /:/, $comps2->[tag_finder('X0', $comps2, 11)];
    	my $x1_pos = tag_finder('X1', $comps2, 11);
    	if($x1_pos) {
	    	(undef, undef, $subopt2) = split /:/, $comps2->[$x1_pos];
    	}
    }
  }
  
  # skip if both are not unique
  if(($xt_1 ne 'U' && $xt_1 ne 'M') && ($xt_2 ne 'U' && $xt_2 ne 'U')) {
  	debug($read_name.' Both reads not SW or unique');
    return;
  }
  
  my ($nm1, $nm2) = (-1, -1); # set as -1 as undef/0 cause extra tests
  if($cigar1 ne '*') { # then 2 must not be
    (undef, undef, $nm1) = split /:/, $comps1->[tag_finder('NM', $comps1, 11)];
  }
  if($cigar2 ne '*') { # then 1 must not be
    (undef, undef, $nm2) = split /:/, $comps2->[tag_finder('NM', $comps2, 11)];
  }
  
  # skip if both are unique with no edit distanct
  # NOTE: added NM check as hides all small indels otherwise
  if($xt_1 eq 'U' && $xt_2 eq 'U' && $nm1 == 0 && $nm2 == 0) {
  	debug($read_name.' Both reads are unique');
    return;
  }
  
  # skip if both have sub optimal hits
  if(($best1 + $subopt1) > 1 && ($best2 + $subopt2) > 1) {
	  debug($read_name.' Both reads suboptimal');
    return;
  }
  
  my $read_length = length $comps1->[9];
  my $max_edit = int ($read_length / $one_edit_per_bps) + 1;
  
  # to be useful for Pindel at least one read must have full match with max of 2 edits
  if($nm1 > $max_edit && $nm2 > $max_edit) {
  	debug($read_name.' Both reads have more than '.$max_edit.' edits (allowing 1 edit per '.$one_edit_per_bps.' bps + 1)');
  	return;
  }
  
  my $uniq;
  
	if($xt_1 && $xt_2) {
		if(($xt_1 eq 'U' || $xt_1 eq 'M') && $nm1 <= $max_edit && $nm2 > 0) {
			$uniq = 1;
		}
		elsif(($xt_2 eq 'U' || $xt_2 eq 'M') && $nm2 <= $max_edit && $nm1 > 0) {
			$uniq = 2;
		}
		else {
			debug($read_name.' Both XT set, nm1 nm2 flags are not of useful values');
			return;
		}
	}
	elsif($cigar1 eq '*' && $nm2 <= $max_edit) {
		$uniq = 2;
	}
	elsif($cigar2 eq '*' && $nm1 <= $max_edit) {
		$uniq = 1;
	}
	else {
		debug($read_name.' One XT, one unmapped, but nm for anchor is >'.$max_edit.' edits (allowing 1 edit per '.$one_edit_per_bps.' bps + 1)');
		return;
	}
	
	my $rec = build_record($comps1, $comps2, $options, $uniq, $read_name);
	my $rec_2;
	if($uniq == 1 && ($xt_2 eq 'U' || $xt_2 eq 'M') && $nm2 <= $max_edit && $nm1 > 0) {
		$uniq = 2; # as previous block will force $uniq=1 if this is true
		$rec_2 = build_record($comps1, $comps2, $options, $uniq, $read_name);
	}
	my @records = ($rec,$rec_2);
  return \@records;
}

sub unmapped_seq {
  my ($unmap_ref, $invert) = @_;
  my $seq = $unmap_ref->[9];
  my $strand = '+';
  if(($unmap_ref->[1] | $strand_mask) == $unmap_ref->[1]) {
    $strand = '-';
  }
  
  if($invert) {
    if($strand eq '+') {
      $seq =~ tr/ACGT/TGCA/;
      $seq = reverse $seq;
    }
  }
  else {
    if($strand eq '-') {
      $seq =~ tr/ACGT/TGCA/;
      $seq = reverse $seq;
    }
  }
  
  if($unmap_ref->[5] eq '*' || $unmap_ref->[5] =~ m/S/) {
  	$seq =~ s/^N+//;
  	$seq =~ s/N+$//;
  }
  
  return $seq;
}

sub length_correct {
  my ($comps) = @_;
  my $length = 0;
  my ($seq, $cigar) = ($comps->[9], $comps->[5]);
  my @cigar_lengths = split /[A-Z]+/, $cigar;
  my @cigar_operations = split /[0-9]+/, $cigar;
  shift @cigar_operations; # will always have empty first element;
  if(@cigar_lengths != @cigar_operations) {
    croak "Cigar string has unmatched lengths and operators for record:\n\t".(join "\t", @{$comps});
  }
  if(@cigar_lengths == 1) {
    if($cigar_lengths[0] == length $seq) {
      $length = $cigar_lengths[0];
    }
    else {
      croak "Cigar string has no variation yet does not match sequence length for:\n\t".(join "\t", @{$comps});
    }
  }
  else {
    for my $i(0..@cigar_lengths-1) {
      if($cigar_operations[$i] eq 'M') {
        $length += $cigar_lengths[$i];
      }
      elsif($cigar_operations[$i] eq 'I') {
        $length += $cigar_lengths[$i];
      }
      elsif($cigar_operations[$i] eq 'S') {
        $length += $cigar_lengths[$i];
      }
      elsif($cigar_operations[$i] eq 'D') {
        $length -= $cigar_lengths[$i];
      }
      else {
        $debug = 1;
        debug('CIGAR: '.$cigar.' seq length: '.length $seq);
        debug('Operation Length:', Dumper(\@cigar_lengths));
        debug('Operation Type:', Dumper(\@cigar_operations));
        croak 'what now? I only understand M,I,S,D for length correction';
      }
    }
  }
  return $length;
}

sub trim_n_from_ends {
	my ($seq) = @_;
	$seq =~ s/^N+//;
	$seq =~ s/N+$//;
	my $seq_len = length $seq; # seq withouth runs of N's at each end
	my $n_in_seq = $seq =~ tr/N//; # returns count of N's within remaning seq
	return ($seq_len, $n_in_seq);
}

sub build_record {
  my ($read1_comps, $read2_comps, $options, $uniq, $read_name) = @_;
  my ($map_comps, $unmap_comps, $pindel_rec);
  
	if($uniq == 1) { # i.e. read 1 is unique
		$map_comps = $read1_comps;
		$unmap_comps = $read2_comps;
	}
	else {
		$map_comps = $read2_comps;
		$unmap_comps = $read1_comps;
	}
  
	my ($seq_len_unmap, $n_in_seq_unmap) = trim_n_from_ends($unmap_comps->[9]);
	my ($seq_len_map, $n_in_seq_map) = trim_n_from_ends($map_comps->[9]);
	
	my $max_n_unmap = int ($seq_len_unmap * $max_n_frac);
	my $max_n_map = int ($seq_len_map * $max_n_frac);
	
	if(	$seq_len_unmap >= $min_seq_hardlimit && $n_in_seq_unmap <= $max_n_unmap
			&&
			$seq_len_map >= $min_seq_hardlimit && $n_in_seq_map <= $max_n_map) {
		my $length = -1; # pindel is 0 based not 1 based so coords of +ve mappings will be -1 to SAM record
		my $strand = '+';
		if(($map_comps->[1] | $strand_mask) == $map_comps->[1]) {
			$strand = '-';
			# -1 so return the location of the last base not the one after
			$length = length_correct($map_comps) - 1;
		}
		
		my $pi_val = $options->{pi} || '*';
		if($options->{'rgs'}) {
			my $rg_pos = tag_finder('RG', $map_comps, 11); # fixed tags are first 11 items so skip
			my (undef, undef, $rg) = split ("/:/xms", $map_comps->[$rg_pos]);
			$pi_val = $options->{'rgs'}->{$rg};
		}
		
		if($map_comps->[2] eq $unmap_comps->[2] && $map_comps->[3] != $unmap_comps->[3] && abs($map_comps->[8]) < $pi_val) {
			# fix for reads that map closely together
			if($strand eq '+') {
				$length -= $pi_val;
			}
			else {
				$length += $pi_val;
			}
		}
		
		$pindel_rec = '@'.$unmap_comps->[0]; # the readname of unmapped record
		if(($unmap_comps->[1] | $first_in_pair_mask) == $map_comps->[1]) {
			$pindel_rec .= "/1";
		}
		elsif(($unmap_comps->[1] | $second_in_pair_mask) == $map_comps->[1]) {
			$pindel_rec .= "/2";
		}
		$pindel_rec .= "\n"; 
		$pindel_rec .= unmapped_seq($unmap_comps)."\n"; # sequence of unmapped read
		
		$pindel_rec .= $strand."\t".$map_comps->[2]."\t"; # refseq name
		
		$pindel_rec .= ($map_comps->[3] + $length)."\t".$map_comps->[4]."\t"; # position, quality
		
		$pindel_rec .= $pi_val."\t"; # predicted insert size
		$pindel_rec .= $options->{'s'}."\n"; # sample tag
		debug($pindel_rec);
	}
	else {
		debug($unmap_comps->[0].' more than allowed Ns in trimmed sequence (allowing max '.($max_n_frac*100).'%) or hard limit of remaning sequence reached - '.$min_seq_hardlimit);
	}
	my @extended_rec;
	if($pindel_rec) {
		@extended_rec = ($pindel_rec, $map_comps->[2], $map_comps->[3]);
	}
  return \@extended_rec;
}

sub tag_finder {
  my ($tag, $comp_ref, $start_at) = @_;
  my $tag_pos;
  for my $i($start_at..@{$comp_ref}-1) {
    my @tag = split ":", $comp_ref->[$i];
    if($tag[0] eq $tag) {
      $tag_pos = $i;
      last;
    }
  }
  return $tag_pos;
}

sub option_builder {
  my ($factory) = @_;

  my %opts = ();

  my $result = GetOptions ( 'i|input=s'   => \$opts{'i'},
			    'o|output=s'  => \$opts{'o'},
			    's|sample=s'  => \$opts{'s'},
			    'c|centromeres=s'  => \$opts{'c'},
			    'om|oldmethod'  => \$opts{'om'},
			    'd|debug'   => \$opts{'d'},
			    'h|help'    => \$opts{'h'},
			    #kb468:
			    'pi=s'      => \$opts{pi},
      );

  usage() if(!$result || !$opts{'i'} || !$opts{'o'} || !$opts{'s'} || $opts{'h'});

  $debug = 1 if($opts{'d'});

  return \%opts;
}

## UTILITIES
sub debug {
  my @data = @_;
  if($debug) {
    foreach my $item(@data) {
      warn $item."\n";
    }
  }
  return;
}

sub usage {

    print <<'USAGE_DOC';

extract_se_map.pl

Description - Extract readpairs where only one end maps fully

Info - It is recommended that you prepare your input BAM file as follows
         samtools view -F 12 -b -o NEW_FILE.bam ORIGINAL.bam
          - remove all reads that are of no use to pindel (both ends unmapped)
         samtools sort -n -m 4000000000 ORIGINAL.bam NEW_namesorted
          - REQUIRED, the above command will need 8GB of RAM

options...
    -i|input        : Input BAM file (req)
    -o|output       : Output filtered ready for pindel
                     (prefix only, _ref[_sw].txt[_a<i>] will be appended, based on anchor read) (req)
    -s|sample       : Sample or label (tumour/normal) (req)
    -d|debug        : Debug messages on and pause at construction of each coord set
    -h|help         : Print this message
    -c|centromeres  : Location of file listing centromere coordinates (gff3 format)
    -om|oldmethod   : Use samtools as converter
    -pi             : Predicted mean if no readgroup is available in the BAM file (kb468)

Outputs - A file with fields relevant for pindel

Example...
  ./bam2pindel_bwa.pl -i /xx/xxxxx.bam -o /yyyyy/yyy -s tumour
  
  For farm use -c can be derived from the array index:
  bsub -o extract.%J.%I -q long -P ??? -J 'extract' "./extract_se_map_comb.pl -i bwa.bam -o bwa -d -s tumour"

Author : kr2

Modified by: kb468

USAGE_DOC

  exit 0;
}
