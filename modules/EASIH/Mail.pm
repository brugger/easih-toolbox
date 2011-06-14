package EASIH::Mail;
# 
# Wrapper for sending std mails using the cam smtp server
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Net::SMTP;


# 
# 
# 
# Kim Brugger (13 Jun 2011)
sub send {
  my ($to, $subject, $content) = @_;

  $to      =~ s/\n//g;
  $subject =~ s/\n//g;

  foreach ( split(/,/, $to) ) {
    
    my $smtp = Net::SMTP->new('ppsw.cam.ac.uk');
    
    my $username = scalar getpwuid $<;
    
    $smtp->mail("$username\@cam.ac.uk");
    
    $smtp->to($_);
    
    $smtp->data();
    $smtp->datasend("Subject: $subject\n");
    $smtp->datasend("To: $_\n");
    $smtp->datasend("\n");
    $smtp->datasend("$content\n");
    $smtp->dataend();
    
    $smtp->quit;
  }
}

1;



