#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

my $project;
my $targ_root;
my @file_type;
my $test;
my $help;

GetOptions ("project=s" => \$project,
	    "ext:s" => \@file_type,
	    "target=s" => \$targ_root,
	    "test" => \$test,
	    "help" => \$help
    );

&help if $help;
my $proj_dir = "/data/$project";
my $targ_dir = "$targ_root/$project";
mkdir $targ_dir unless -e $targ_dir;

my %file2loctn = ('fq' => 'SEQUENCE',
		  'fq.gz' => 'SEQUENCE',
		  'bam'=> 'ALIGNMENT',
		  'vcf'=> 'VARIATION',
		  'csv'=> 'VARIATION',
		  'pdf'=> 'QC',
		  'flagstat'=> 'QC',
		  'bai'=> 'ALIGNMENT'
);
my %file2src = ('fq' => 'raw',
		'fq.gz' => 'raw',
		'bam'=> 'MARIS',
		'vcf'=> 'MARIS',
		'csv'=> 'MARIS',
		'pdf'=> 'raw',
		'flagstat'=> 'MARIS',
		'bai'=> 'MARIS'
    );

if($file_type[0] eq 'all') {
    @file_type = keys %file2loctn;
}

foreach (@file_type) {
    if($file2src{$_}) {
	my $search = "$proj_dir/".$file2src{$_}."/*".$_;
	my @list = glob("$search");
	foreach my $file (@list) {
	    my($filename, $directories) = fileparse($file);
	    my $dest = "$targ_dir/".$file2loctn{$_};
	    mkdir $dest unless (-e $dest);
	    defined $test ? system("echo cp $file \"$dest/$filename\"") : system("cp $file \"$dest/$filename\"");
	}
    }
    else {
	warn "$_ is invalid filetype : NOT COPIED\n";
    }
}

sub help {
    



}

=pod

=head1 prepare_data.pl

=over 

=item Script to copy results files to standard directory structure for return to customer.

=back

=head2 Usage

=over

=item -project  EASIH project code eg A01

=item -target  root dir to copy files eg /media/USBdrive

=item -ext  file extensions to copy eg -ext fq.gz -ext pdf

=item -help show this help

=item -test  just output commands to STDOUT rather than execute them

=back 

=head2 Example

=over 

=item prepare_data.pl -project Z01 -target /media/KINGSTON -ext bam -ext vcf

=item will copy bam (and bai) and vcf files to a directories under a project folder called Z01 on the mounted KINGSTON USB drive

=back

=cut
