package EASIH::HTML::style;
#!/usr/bin/perl 
# 
# Based my code from the mutagen genome browser 
# 
# 
# Kim Brugger (15 Jul 2011), contact: kim.brugger@easih.ac.uk

use strict;

package html::style;
require Exporter;
require AutoLoader;

# set the version for version checking
my $VERSION     = "0.01";

my @ISA = qw(Exporter AutoLoader);

our @EXPORT = qw(check_user);

@ISA = qw(Exporter AutoLoader);

@EXPORT = qw();

# headers and simple text styles.

sub h1 {
  my ($text) = @_;
  return "<H1>$text</H1>";
}

sub h2 {
  my ($text) = @_;
  return "<H2>$text</H2>";
}

sub h3 {
  my ($text) = @_;
  return "<H3>$text</H3>";
}

sub h4 {
  my ($text) = @_;
  return "<H4>$text</H4>";
}

sub h5 {
  my ($text) = @_;
  return "<H5>$text</H5>";
}

sub h6 {
  my ($text) = @_;
  return "<H6>$text</H6>";
}

sub bold {
  my ($text) = @_;
  return "<B>$text</B>";
}

sub italics {
  my ($text) = @_;
  return "<I>$text</I>";
}

sub monospace {
  my ($text) = @_;
  return "<TT>$text</TT>";
}

sub underline {
  my ($text) = @_;
  return "<U>$text</U>";
}

sub hr {
  my ($width) = @_;
  return "<HR width='$width'>\n" if $width;
  return "<HR>\n";
}

sub center {
  my ($text) = @_;
  return "<CENTER> $text </CENTER>\n";
}

sub break {
  return "<BR>\n";
}

1;
