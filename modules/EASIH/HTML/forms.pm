package EASIH::HTML::forms;
#!/usr/bin/perl 
# 
# Based my code from the mutagen genome browser 
# 
# 
# Kim Brugger (15 Jul 2011), contact: kim.brugger@easih.ac.uk

use strict;

package html::forms;
require Exporter;
require AutoLoader;

# set the version for version checking
my $VERSION     = "0.01";

my @ISA = qw(Exporter AutoLoader);

our @EXPORT = qw(check_user);

@ISA = qw(Exporter AutoLoader);

@EXPORT = qw();

sub text {
  my ($name, $value, $readonly, $size, $maxlength) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  
  $call_hash{'name'}      = $name;
  $call_hash{'type'}      = "text";
  $call_hash{'value'}     = $value     if ($value);
  $call_hash{'readonly'}  = $readonly  if ($readonly);
  $call_hash{'size'}      = $size      if ($size);
  $call_hash{'maxlength'} = $maxlength if ($maxlength);

  return(html::generic_form_element(\%call_hash));
}

sub password {
  my ($name, $value, $size, $maxlength) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  
  $call_hash{'name'}      = $name;
  $call_hash{'type'}      = "password";
  $call_hash{'value'}     = $value     if ($value);
  $call_hash{'size'}      = $size      if ($size);
  $call_hash{'maxlength'} = $maxlength if ($maxlength);

  return(html::generic_form_element(\%call_hash));
};

sub checkbox {
  my ($name, $values, $default) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  return 0 if (ref($values) ne "ARRAY");

  my $return_string = "";
  foreach my $value (@$values) {

    $call_hash{'name'}      = $name;
    $call_hash{'type'}      = "checkbox";
    $call_hash{'value'}     = $value    if ($value);
    $call_hash{'default'}   = $default   if ($default);

    $return_string .= html::generic_form_element(\%call_hash);
  }

  return($return_string);
};

sub checkbox_group {
  my ($name, $values, $default) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  return 0 if (ref($values) ne "ARRAY");

  my $return_string = "";
  foreach my $value (@$values) {

    $call_hash{'name'}      = $name;
    $call_hash{'type'}      = "checkbox";
    $call_hash{'value'}     = $value    if ($value);
    $call_hash{'default'}   = $default   if ($default);

    $return_string .= html::generic_form_element(\%call_hash);
  }

  return($return_string);
};

sub radio {
  my ($name, $values, $default) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  return 0 if (ref($values) ne "ARRAY");

  my $return_string = "";
  foreach my $value (@$values) {

    $call_hash{'name'}      = $name;
    $call_hash{'type'}      = "radio";
    $call_hash{'value'}     = $value    if ($value);
    $call_hash{'default'}   = $default   if ($default);

    $return_string .= html::generic_form_element(\%call_hash);
  }

  return($return_string);
};

sub submit {
  my ($name, $value) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  
  $call_hash{'name'}      = $name;
  $call_hash{'type'}      = "submit";
  $call_hash{'value'}     = $value     if ($value);

  return(html::generic_form_element(\%call_hash));
};

sub reset {
  my ($name, $value) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  
  $call_hash{'name'}      = $name;
  $call_hash{'type'}      = "reset";
  $call_hash{'value'}     = $value     if ($value);

  return(html::generic_form_element(\%call_hash));
};


sub hidden {
  my ($name, $value) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  
  $call_hash{'name'}      = $name;
  $call_hash{'type'}      = "hidden";
  $call_hash{'value'}     = $value     if ($value);

  return(html::generic_form_element(\%call_hash));
};

sub button{
  my ($name, $value) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  
  $call_hash{'name'}      = $name;
  $call_hash{'type'}      = "button";
  $call_hash{'value'}     = $value     if ($value);

  return(html::generic_form_element(\%call_hash));
};

sub textarea {
  my ($name, $value, $readonly, $rows, $cols, $maxlength) = @_;

  my %call_hash = ();

  return 0 if (! $name);
  
  $call_hash{'name'}      = $name;
  $call_hash{'type'}      = "textarea";
  $call_hash{'value'}     = $value     if ($value);
  $call_hash{'readonly'}  = $readonly  if ($readonly);
  $call_hash{'rows'}      = $rows      if ($rows);
  $call_hash{'cols'}      = $cols      if ($cols);
  $call_hash{'maxlength'} = $maxlength if ($maxlength);

  return(html::generic_form_element(\%call_hash));
};


END { 

} 

1;
