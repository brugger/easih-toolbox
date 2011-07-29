package EASIH::HTML;
#!/usr/bin/perl 
# 
# Based my code from the mutagen genome browser 
# 
# 
# Kim Brugger (15 Jul 2011), contact: kim.brugger@easih.ac.uk

use strict;

use EASIH::HTML::forms;
use EASIH::HTML::style;

use CGI;

require Exporter;
require AutoLoader;

# set the version for version checking
my $VERSION     = "0.01";

my @ISA = qw(VERSION Exporter AutoLoader);

our @EXPORT = qw();
@ISA = qw(Exporter AutoLoader);
@EXPORT = qw();

#exported variables, this will be sweet.
use vars qw (%parameters $query $remote_host);


BEGIN {
  $query = new CGI;
  $remote_host = $query->remote_host();
  %parameters = $query->Vars;

  foreach my $key (keys %parameters) {
    my @dbs = split("\0",$parameters{$key});
    if (@dbs > 1) {
      $parameters{$key} = \@dbs;
    }
  }
}

sub head {
  my $query = new CGI;
  return $query->header;
  
}

sub start {
  my ($title, $css, $values, $no_start) = @_;
  
  my $return_string = "";
  $return_string = head() if (!$no_start);
  $return_string .= "<HTML>\n";
  $return_string .= "<TITLE>$title</TITLE>\n" if ($title);
  $return_string .= "<LINK rel='stylesheet' href='$css' type='text/css' />\n" if ($css);

  $return_string .= "<BODY";

  foreach my $key (keys %$values) {
    $return_string .= " $key='$$values{$key}'";
  }
  return "$return_string >\n"; 
}

sub end {
  return  "</html>";
}


sub link {
  my ($link, $text, $target) = @_;

  return "<A HREF='$link' TARGET='$target' >$text</A>" if $target;
  return "<A HREF='$link' >$text</A>";
}

sub image {
  
}

sub redirect {
  my ($url, $time) = @_;
  

  return "<HEAD><META HTTP-EQUIV='Refresh'  CONTENT='$time'; URL='$url'></HEAD>";
}

# Creates a form element from a hash table, this will be wrapped by other
# functions but it is nice to have a generic function for this.
# 
# Kim Brugger (20 Oct 2003)
sub generic_form_element {
  my ($element) = @_;
  
  # Following fields and values are accepted in the element-hash:
  # 
  # type:      text|password|checkbox|radio|submit|reset|file|hidden|image|button|label|menu|popup
  # name:      field name
  # value:     value for the element.
  # values:    for a menu or popup this is several values in an array.
  # defaults:  what has been selected (if any) (an array)
  # default:   what has been selected (if any) (single value)
  # multiple:  possible to select several fields or not
  # readonly:  only for text fields
  # checked:   only for radio and check boxes
  # rows:      controls the size of a textarea
  # cols:      controls the size of a textarea
  # size:      the size of the element
  # maxlength: man length of the input element
  # 
  # Additional (just for this implementation of the HTML form:
  #  pre_string: added in front of the INPUT element
  #  post_string: added in after the INPUT element
  
  my $return_string = "";;
  $return_string .= "$$element{pre_string}" if $$element{pre_string};

  $return_string .= "<INPUT type='$$element{type}'" if ($$element{type} ne "textarea" &&
							$$element{type} ne "menu" &&
							$$element{type} ne "popup");
  $return_string .= "<TEXTAREA" if ($$element{type} eq "textarea");
  $return_string .= "<SELECT" if ($$element{type} eq "menu" || $$element{type} eq "popup");

  $$element{size} = 1 if ($$element{type} eq "popup");

  foreach my $key (keys %$element) {
    next if ($key eq 'type');
    next if ($key eq 'pre_string');
    next if ($key eq 'post_string');
    next if ($key eq 'values');
    next if ($key eq 'labels');
    next if ($key eq 'defaults');
    next if ($key eq 'default');
    next if ($key eq 'value' && $$element{type} eq "textarea");
    next if (not defined $key);

    $return_string .= " $key='$$element{$key}'" if ($$element{$key});
  }	

  $return_string .= ">\n";

  # if an text area the value have to be placed between the textarea tags.
  if ($$element{type} eq "textarea") {
    $return_string .= "$$element{value}\n" if ($$element{value});
    $return_string .= "</TEXTAREA>\n";
  }

  # Menu's and popups are even worse, so here we go...
  if ($$element{type} eq "menu" || $$element{type} eq "popup") {
    # first the default selections are placed in a hash.
    my %defaults = ();
    foreach my $value (@{$$element{defaults}}) {
      next if (not defined $value);
      $defaults{$value} = 1;

    }
    $defaults{$$element{'default'}} = 1 if ($$element{default});

    # then the tags are constructed.
    foreach my $value (@{$$element{values}}) {
      $return_string .= "  <OPTION ";
      
      $return_string .= "selected " if ($defaults{$value});
      $return_string .= "value='$value'>";

      $return_string .= "$$element{labels}{$value}</OPTION>\n" if ($$element{labels}{$value});
      $return_string .= "$value</OPTION>\n" if (!$$element{labels}{$value});
    }
    $return_string .= "</SELECT>\n";
  }


  $return_string .= " $$element{post_string}" if $$element{post_string};

#  print STDERR  "$return_string\n";
  return "$return_string\n";
}

sub make_form {
  my ($elements, $url, $method) = @_;

  my $return_string = start_form($url, $method);

  foreach my $element (@$elements) {
    $return_string .= generic_form_element($element);
  }

  $return_string .= end_form();
  
#  print STDERR "$return_string\n\n\n";

  return $return_string;
}

sub start_form {
  my ($url, $method, $target) = @_;
  
  my $return_string = "<FORM";
  $return_string .= " action='$url'" if ($url);
  $return_string .= " method='$method'" if ($method);
  $return_string .= " method='POST'" if (!$method);
  $return_string .= " target='$target'" if ($target);
  

  return "$return_string>\n";
}


sub end_form {
  return "</FORM>\n"
}


sub start_multiform {
  my ($url, $method, $target) = @_;
  
  my $return_string = "<FORM enctype='multipart/form-data'";
  $return_string .= " action='$url'" if ($url);
  $return_string .= " method='$method'" if ($method);
  $return_string .= " method='POST'" if (!$method);
  $return_string .= " target='$target'" if ($target);
  

  return "$return_string>\n";
}


sub save_upload_file {
  my ($identifier) = @_;

  use POSIX qw( tmpnam );

  my $tmpout =  tmpnam();

  open FIL, ">$tmpout" || die "Could not open '$tmpout': $!\n";

  # return list of errors if any
  
  my ( $fileobj, $fh, $type, $cd, $content );
  
  $fileobj = $query->param($identifier);
  
  if ( $fileobj ) {
    $fh = $query->upload($identifier);
    
    $type = $query->uploadInfo( $fileobj )->{ 'Content-Type' };

    if ( $type =~ m!text/plain! || $type =~ m!pubmed/text! || $type =~ m!application/octet-stream!) {
      $cd       = $query->uploadInfo( $fileobj )->{ 'Content-Disposition' };

      undef $/;
      $content = <$fh>;
      print FIL $content;
      $/ = "\n";
      
      if ( $query->cgi_error ) {
	print STDERR $query->cgi_error;
	return undef;
      }
      
    }
    else {
      print STDERR "Upload file type must be 'text/plain' not '$type'\n";
      return undef;
    }
  }

  close FIL || die "Could not close '$tmpout': $!\n";
  
  return $tmpout;
}


#
# Creates an advanced table, the function expects an array of arrays, and a border or not flag.
# 
# Kim Brugger (20 Oct 2003)
sub advanced_table {
  my ($cells,      # the cells as an array of arrays of hashes, named according to the language:
                   #              %cell = (value, bgcolor, width, colspan, 
                   #                       rowspan, cellhalign, cellvalign);
      $border,     # border or not.
      $padding,    # how big the cells are (padded around the text).
      $spacing,    # how the cells should be spaced
      $bgcolour,   # the colour of the table.
      $spanning,   # Wether there are spanning cells, or if they should be "padded".
      $tablewidth, # how wide the table should be, this is a string, 
                   # so we can handle both pixel width and percentages
      $class       # Class so CSS can define the behaviour
      ) = @_;

  my $width = 0;
  foreach my $row (@$cells) {
    $width = @$row if ($width < @$row);
  }
  
  my $return_string = table_start($border, $padding, $spacing, $bgcolour, $tablewidth, $class);
  
  foreach my $row (@$cells) {
    
    $return_string .= "  <TR>";
    for (my $i=0; $i<$width;$i++) {
      
      if ($$row[$i]) {
	$return_string .= "<TD";
	
	if (ref $$row[$i] eq "HASH") {

	  foreach my $key (keys %{$$row[$i]}) {
	    next if ($key eq 'value');
	    $return_string .= " $key='$$row[$i]{$key}'";
	  }	
	  $return_string .= ">$$row[$i]{'value'}</TD>" if ($$row[$i]{'value'});
	}
	else {
	  $return_string .= ">$$row[$i]</TD>" if ($$row[$i]);
	}
      }
      else {
	$return_string .= "<TD>&nbsp;</TD>" if (!$spanning);
      }

    }
    $return_string .= "</TR>\n";
  }

  $return_string .= table_end();

  return $return_string;
}


#
# Creates a simple table, the function expects an array of arrays, and a border or not flag.
# 
# Kim Brugger (20 Oct 2003)
sub table {
  my ($cells,     # the cells as an array of arrays of values.
      $border,    # border or not.
      $padding,   # how big the cells are (padded around the text).
      $spacing,   # how the cells should be spaced
      $bgcolour,  # the colour of the table.
      $tablewidth, # how wide the table should be, this is a string, so we can handle both pixel width and percentages
      ) = @_;


  my $width = 0;
  foreach my $row (@$cells) {
    $width = @$row if ($width < @$row);
  }

  my $return_string = table_start($border, $padding, $spacing, $bgcolour, $tablewidth);

  foreach my $row (@$cells) {
    
    $return_string .= "  <TR>";
    for (my $i=0; $i<$width;$i++) {
      $$row[$i] = "&nbsp;" if (not defined $$row[$i]);
      $return_string .= "<TD>$$row[$i]</TD>" 
    }
    $return_string .= "</TR>\n";
  }

  $return_string .= table_end();

  return $return_string;
}


sub table_start {
  my ($border, $padding, $spacing, $bgcolour, $width, $class) = @_;
  
  my $return_string .= "<TABLE";
  $return_string .= " border='$border'"        if $border;
  $return_string .= " cellspacing='$spacing'"  if $spacing;
  $return_string .= " cellpadding='$padding'"  if $padding;
  $return_string .= " bgcolor='$bgcolour'"     if $bgcolour;
  $return_string .= " width='$width'"          if $width;
  $return_string .= " class='$class'"          if $class;
  $return_string .= " >\n";

  return $return_string;
}

sub table_end {
  return "</TABLE>\n";

}

sub list {
  my ($values) = @_;
  return "" if (ref($values) ne "ARRAY");

  my $return_string = "<UL>\n";
  foreach my $value (@$values) {
    $return_string .= "<LI> $value\n";
  }
 $return_string .= "</UL>\n";

  return $return_string;
}

sub ordered_list {
  my ($values) = @_;
  return "" if (ref($values) ne "ARRAY");

  my $return_string = "<DL>\n";
  foreach my $value (@$values) {
    $return_string .= "<LI> $value\n";
  }
 $return_string .= "</DL>\n";

  return $return_string;
}

sub contact {
  my ($name, $email) = @_;

  return "<CENTER><A href='mailto:$email'>$name</A></CENTER>";
}

sub checkbox_table {
  my ($element, $width) = @_;
  

  my $labels = $$element{'labels'};

  my $html = "";
  my @rows = ();
  my @cells= ();
  $width = 4 if (!$width);
  for(my $i=0;$i < @{$$element{values}}; $i++) {
    my %call_hash = (type => "checkbox",
		     name => $$element{'name'},
		     defaults=> $$element{'defaults'},
		     value=>$$element{values}[$i]
		     );
    push @rows, generic_form_element(\%call_hash) . "$$labels{$$element{values}[$i]}";
    
    if (($i+1)%$width == 0) {
      my @rs = @rows;
      push @cells, \@rs;
      @rows = ();
    }
  }

  push @cells, \@rows if (@rows);

  return table(\@cells, 0, 1, 1);
}

sub dump_params {
  use Data::Dumper;
  print STDERR "---- EASIH::HTML::parameters ----\n";
  print STDERR Dumper (\%parameters);
  print STDERR "--------------------------\n";
} 

1;
