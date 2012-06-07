#!/usr/bin/perl 

# Create a set of multiplexed sample tags which exhibit specified
# error correction and detection properties.

# Tags of specified length are generated randomly, then either
# accepted or rejected. If a newly generated tag collides with a
# previously accepted tag in a way which would invalidate the error
# correction guarantees, it is rejected. There is also provision for
# applying arbitrary criteria based on the tag itself, such as base
# content.

# Two tables are maintained. The correction_table maps all
# perturbations of the accepted tags back to the correct tags. This
# would be used an any error-correction algorithnm. The
# detectable_errors table contains the larger "shadow" of the accepted
# tags: those tags which must be rejected because their correction
# table would include a tag which is already a "detectable error" case
# for an previously accepted tag.

# A newly generated tag is rejected if its correction table overlaps
# with the cumulative detectable_errors table of all previously
# accepted tags.

# Tom Skelly     ts6@sanger.ac.uk

use strict;
use warnings;
use English qw(-no_match_vars);
use Carp;
use Getopt::Long;

my $VERSION = '20100216.01';

sub process;
sub usage;
sub initialise;
sub generate_tag($);
sub tag_by_index($$);
sub validate_tag($);
sub validate_base_content($);
sub generate_collision_set($$);
sub intersect($$);
sub add_to_table($$);
sub off_by_one($);
sub off_by_two($);
sub off_by_shift($);
sub validate_edit_distance($);
sub validate_shift_distance($);

my $DEFAULT_LEN = 7;
my $DEFAULT_NUM = 512;
my $DEFAULT_WINDOW = 100;

my @BASE_TAGS = qw[A C G T];

my $opts = initialise;
my $verbose = $opts->{verbose} || 0;

process;   # execute the main routine

exit;

# ----------------------------------------------------------------------
sub usage {

  ## no critic

  my ($short_name) = $PROGRAM_NAME =~ /([^\/]+)$/xm;

  print STDERR "\n";
  print STDERR "$short_name version $VERSION\n";
  print STDERR "\n";
  print STDERR "    options:\n";
  print STDERR "\n";
  print STDERR "    --help         print this message and quit\n";
  print STDERR "    --verbose      write step-by-step stats to STDERR\n";
  print STDERR "    --num          number of tags to (try to) generate\n";
  print STDERR "    --len          length of tag\n";
  print STDERR "    --file         read candidate tags from file, rather than randomly generating them\n";
  print STDERR "    --add-tags     add random tags to set after reading file\n";
  print STDERR "    --all-tags     try all possible (rather than random) tags of specified length\n";
  print STDERR "    --single       correct/detect single-base errors\n";
  print STDERR "    --double       correct/detect two-base errors\n";
  print STDERR "    --shift        correct/detect single-shift errors\n";
  print STDERR "    --min-bases    minimum count of each base in tag\n";
  print STDERR "    --seed         random seed for tag generation\n";
  print STDERR "    --window       quit early if no hits in this many tries (def $DEFAULT_WINDOW)\n";
  print STDERR "\n";
  print STDERR "  tag set is output to STDOUT\n";
  print STDERR "  single, double, and shift options accept 'correct' or 'detect'\n";
  print STDERR "  if no seed is specified, a random seed will be chosen\n";
  print STDERR "\n";
  print STDERR "example:\n";
  print STDERR "\n";
  print STDERR "  $short_name --verbose --len 7 --single correct --double detect --shift detect\n";
  print STDERR "\n";

}

# ----------------------------------------------------------------------
sub initialise {

  my %opts;
  my $rc = GetOptions(\%opts,
                      'help',
                      'verbose',
                      'num:i',
                      'len:i',
                      'file:s',
                      'add-tags',
                      'all-tags',
                      'single:s',
                      'double:s',
                      'shift:s',
                      'min-bases:i',
                      'seed:i',
                      'window:i',
                      );

  if ( ! $rc) {
    print {*STDERR} "\nerror in command line parameters\n" or croak 'print failed';
    usage;
    exit 1;
  }

  if (exists $opts{'help'}) {
    usage;
    exit;
  }

  foreach my $parm (qw[single double shift]) {
    if (exists $opts{$parm} and $opts{$parm} ne 'detect' and $opts{$parm} ne 'correct') {
      print STDERR "\n\nERROR: invalid value specified for $parm paramater\n";
      usage;
      exit;
    }
  }

  $opts{'seed'} ||= rand 999999;    # want to know this explicitly, so we can print it later
  srand $opts{'seed'};

  return \%opts;

}

# ----------------------------------------------------------------------
sub process {

  my @tags;                     # accepted tags
  my $correction_table  = {};   # correction table for accepted tags (see note below)
  my $detectable_errors = {};   # table of detectable errors

  # counters for determining when to quit

  my $accept_count = 0;
  my $window_count = 0;
  my $accept_in_window = 0;

  my $num_tags = $opts->{'num'}      || $DEFAULT_NUM;
  my $window   = $opts->{'window'}   || $DEFAULT_WINDOW;

  my $new_tag;
  my $tag_len;

  my $fh;
  if (exists $opts->{'file'}) {
    open $fh, $opts->{'file'} or croak "could not open input file ".$opts->{'file'};
    $tag_len = $opts->{'len'} || 0;               # let tags in file determine length if not specified
  } else {
    $tag_len = $opts->{'len'} || $DEFAULT_LEN;    # use default length if not specified
  }

  my $all_tags = $opts->{'all-tags'} || 0;
  my $last_tag = 4 ** $tag_len;
  my $next_tag = 0;

  while ($accept_count < $num_tags) {

    if ($fh) {

      $new_tag = read_tag ($fh, $tag_len);                            # get candidate tag from file

      if ( ! $new_tag) {                                              # end of file?
        $verbose && printf STDERR "\nread %d tags from file %s\n", $., $opts->{'file'};
        close $fh;
        $fh = '';
        last unless exists $opts->{'add-tags'};                       # quit unless asked to add random tags
        $accept_in_window = 0;                                        # restart the window counter
        $new_tag = generate_tag ($tag_len);
        $verbose && print STDERR "continuing with randomly generated tags\n\n";
      }

    } elsif ($all_tags) {
      $new_tag = tag_by_index ($next_tag, $tag_len);
      last if ++$next_tag >= $last_tag;
    } else {
      $new_tag = generate_tag ($tag_len);                             # randomly generate a candidate tag
    }

    my $corrections = generate_collision_set ($new_tag, 'correct');   # errors we want to correct

    # If the candidate tag's correction set has no overlap with the
    # set of correctable+detectable errors from all previously
    # accepted tags, and also passes the ad-hoc checks, accept it.

    if ( ! intersect ($corrections, $detectable_errors)) {
      
      my $known_tag = grep(/$new_tag/, @tags);
      
      if (! $known_tag && validate_tag ($new_tag)) {

        push @tags, $new_tag;

        # %correction_table turns out not to be actually used for
        # anything in this processing. We could remove it, and all
        # references to it. It would be needed in the correction
        # algorithm -- but it could be re-generated from the tag set
        # and error-correction spec. I've only kept it on here in case
        # we want to print it at some point.

        add_to_table ($correction_table, {$new_tag => "$new_tag:tag"});
        add_to_table ($correction_table, $corrections);

        my $detectables = generate_collision_set ($new_tag, 'detect');

        add_to_table ($detectable_errors, {$new_tag => "$new_tag:tag"});
        add_to_table ($detectable_errors, $corrections);
        add_to_table ($detectable_errors, $detectables);

        ++$accept_count;
        ++$accept_in_window;

      }

    }

    if ( ! $fh and ! $all_tags) {                           # don't quit if we're reading a file
      if (++$window_count >= $window) {                     # check whether we're still finding tags
        $verbose && printf STDERR "%3d of %3d accepted\n", $accept_in_window, $window;
        last if $accept_in_window == 0;
        $window_count = 0;
        $accept_in_window = 0;
      }
    }

  }

  validate_edit_distance (\@tags);
  validate_shift_distance (\@tags);

  print join ("\n", @tags), "\n";

  if ( ! exists $opts->{'file'} or exists $opts->{'add-tags'}) {
    printf STDERR "\nfound %d tags with random seed %d\n\n", scalar @tags, $opts->{'seed'};
  } else {
    printf STDERR "\nfound %d tags\n\n", scalar @tags;     # didn't use seed
  }

}

# ----------------------------------------------------------------------
# Randomly generate a tag of the required length.

sub generate_tag ($) {

  my ($tag_len) = @_;
  my @bases;
  
  my $check_sum = 0;

  foreach my $i (1..$tag_len) {
    my $rand_ix = int(rand(@BASE_TAGS));
    $check_sum += $rand_ix;

    push @bases, $BASE_TAGS[$rand_ix];
  }

  $check_sum %= 4;
  push @bases, $BASE_TAGS[$check_sum];
  

  return join '', @bases;

}

# ----------------------------------------------------------------------
# Generate a tag given an index into the complete set of tags of a
# given length.

sub tag_by_index ($$) {

  my ($index, $tag_len) = @_;
  my @bases;

  foreach my $i (1..$tag_len) {
    my $ix = $index & 0x0003;
    $index >>= 2;
    unshift @bases, $BASE_TAGS[$ix];
  }

  return join '', @bases;

}

# ----------------------------------------------------------------------
# Read the next candidate tag from a file.

sub read_tag ($) {

  my ($fh, $tag_len) = @_;

  my $tag = <$fh>;

  return '' unless $tag;

  chomp $tag;
  if ($tag_len > 0) {
    $tag = substr $tag, 0, $tag_len;
  }

  return uc $tag;    # uppercase it

}

# ----------------------------------------------------------------------
# Validate a tag against an arbitrary set of criteria which depend
# only on the tag itself, not on the prior history of accepted tags --
# e.g., distribution of bases in the tag.

sub validate_tag ($) {

  my ($tag) = @_;
  my $ret = 1;

  $ret &&= validate_base_content ($tag);    # currently this is the only check

  return $ret;

}

# ----------------------------------------------------------------------
# Verify that a tag contains at least N of each of the 4 bases.

sub validate_base_content ($) {

  my ($tag) = @_;
  my $ret = 1;

  return $ret unless exists $opts->{'min-bases'};

  my $min = $opts->{'min-bases'};
  my %counts = map {$_ => 0} @BASE_TAGS;

  foreach my $base (split '', $tag) {
    ++$counts{$base};
  }

  for my $base (@BASE_TAGS) {
    if ($counts{$base} < $min) {
      $ret = 0;
      $verbose && print STDERR "$tag fails min-bases\n";
      last;
    }
  }

  return $ret;

}

# ----------------------------------------------------------------------
# Given a tag, generate the list of all colliding tags, according to
# some set of criteria as specified in the command line
# options. Currently suported options include single, double and shift
# errors. Return a ref to a hash keyed by the colliding tag. Hash
# value is a string containing the original tag and the collision
# criterion, used only for debugging and validation.

sub generate_collision_set ($$) {

  my ($tag, $criterion) = @_;
  my $new_set = {};

  if (exists $opts->{'single'} and $opts->{'single'} eq $criterion) {
    foreach my $new ( @{ off_by_one($tag) } ) {
      $new_set ->{$new} = "$tag:off-1";
    }
  }

  if (exists $opts->{'double'} and $opts->{'double'} eq $criterion) {
    foreach my $new ( @{ off_by_two($tag) } ) {
      $new_set ->{$new} = "$tag:off-2";
    }
  }

  if (exists $opts->{'shift'} and $opts->{'shift'} eq $criterion) {
    foreach my $new ( @{ off_by_shift($tag) } ) {
      $new_set ->{$new} = "$tag:shift";
    }
  }

  return $new_set;

}

# ----------------------------------------------------------------------
# Given two hash refs, return a boolean indicating whether they have
# any keys in common. The keys of the first hash will be searched for
# in the second hash, so specifying the smaller hash as the first
# argument will be faster.

sub intersect ($$) {

  my ($hash_1, $hash_2) = @_;
  my $ret = 0;

  foreach my $key (keys %$hash_1) {
    if (exists $hash_2->{$key}) {
      $verbose && print STDERR $hash_2->{$key}, ' -> ', $key, "\n";
      $ret = 1;
      last;
    }
  }

  return $ret;

}

# ----------------------------------------------------------------------
# Given two hash refs, add the entries of the second hash to the first
# hash. Existing entries in the first hash will be overwritten by any
# matching entries from the second hash. Return the first hash ref,
# which now points to the updated hash.

sub add_to_table ($$) {

  my ($target, $new_ents) = @_;

  foreach my $new (keys %$new_ents) {
    $target->{$new} = $new_ents->{$new};
  }

  return $target;

}

# ----------------------------------------------------------------------
# Given a word, return the list of all words which differ from it in
# one base position.

sub off_by_one ($) {

  my ($word) = @_;
  my @wrong_words;

  for (my $i=0; $i<length($word); ++$i) {

    for my $error (@BASE_TAGS) {

      next if substr ($word, $i, 1) eq $error;

      my $new_word = $word;
      substr ($new_word, $i, 1) = $error;
      push @wrong_words, $new_word;

    }
  }

  return \@wrong_words;

}

# ----------------------------------------------------------------------
# Given a word, return the list of all words which differ from it in
# two base positions.

sub off_by_two ($) {

  my ($word) = @_;
  my @wrong_words;

  for (my $i=0; $i<length($word); ++$i) {

    for my $error (@BASE_TAGS) {

      next if substr ($word, $i, 1) eq $error;

      my $new_word = $word;
      substr ($new_word, $i, 1) = $error;

      for (my $j=$i+1; $j<length($word); ++$j) {

        for my $error_2 (@BASE_TAGS) {

          next if substr ($new_word, $j, 1) eq $error_2;

          my $newer_word = $new_word;
          substr ($newer_word, $j, 1) = $error_2;
          push @wrong_words, $newer_word;

        }

      }
    }
  }

  return \@wrong_words;

}

# ----------------------------------------------------------------------
# Given a word, return the list of all words which differ from it by a
# one-base shift in either direction, with all 4 possible bases shifted in.

sub off_by_shift ($) {

  my ($word) = @_;
  my @wrong_words;
  my %prior;

  my @bases = split //, $word;

  for my $error (@BASE_TAGS) {

    my $new_left  = join '', @bases[1..@bases-1], $error;
    my $new_right = join '', $error, @bases[0..@bases-2];

    if ( ! exists $prior{$new_left}) {
      ++$prior{$new_left};
      push @wrong_words, $new_left;
    } else {
      $verbose && print STDERR "duplicate $word -> $new_left\n";
    }

    if ( ! exists $prior{$new_right}) {
      ++$prior{$new_right};
      push @wrong_words, $new_right;
    } else {
      $verbose && print STDERR "duplicate $word -> $new_right\n";
    }

  }

  return \@wrong_words;

}

# ----------------------------------------------------------------------
# For each codeword, print a histogram of edit (Hamming) distances
# from all other codewords. Then following table lists the minimum
# distance required to support various error-correction options that
# can be specified on the command line:

#     2  single-error detect
#     3  single-error correct
#     4  single-error correct, double error detect
#     5  single-error correct, double error correct

sub validate_edit_distance ($) {

  my ($codes) = @_;

  print STDERR "\nValidating minimum edit distance\n\n";

  foreach my $word (@$codes) {

    my @histo = (0) x (length($word)+1);

    foreach my $other_word (@$codes) {

      next if $word eq $other_word;

      my $misses = 0;

      for (my $i=0; $i<length($word); ++$i) {
        if (substr($word, $i, 1) ne substr($other_word, $i, 1)) {
          ++$misses;
        }
      }

      ++$histo[$misses];

    }

    print STDERR "$word @histo\n";

  }

}

# ----------------------------------------------------------------------
# Verify that no pair of codewords differ by only a shift. This will
# be true only if the --shift command line paramater was specified.

sub validate_shift_distance ($) {

  my ($codes) = @_;
  my $hits = 0;

  print STDERR "\nValidating shift distance\n\n";

  my $len = length $codes->[0];

  foreach my $ix1 (0..@$codes-2) {

    my $word1 = $codes->[$ix1];

    foreach my $ix2 ($ix1+1..@$codes-1) {

      my $word2 = $codes->[$ix2];

      if (substr ($word1, 0, $len-1) eq substr ($word2, 1, $len-1)) {
        print STDERR "$word1 << $word2\n";
        ++$hits;
      }
      if (substr ($word1, 1, $len-1) eq substr ($word2, 0, $len-1)) {
        print STDERR "$word1 >> $word2\n";
        ++$hits;
      }

    }

  }

  printf STDERR "\n%d hits found\n\n", $hits;

}
