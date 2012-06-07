#!/usr/bin/perl -w
#*************************************************************************
#
#   File:       DataUploadToHPC.pl
#   
#   Version:    V1.2
#   Date:       07.06.11
#   Function:   Uploading sequencing data to HPC
#   
#   Copyright:  (c) University of Cambridge / S.V.V. Deevi 2011
#   Author:     S.V.V. Deevi
#   Address:    Eastern Sequence Informatics Hub (EASIH),
#               The University of Cambridge,
#               Level 4, Laboratory Block Box 268,
#               Addenbrooke's Hospital,
#               Cambridge,
#               CB2 0QQ,
#               UK.
#   EMail:      svvd2@cam.ac.uk (or) sri.deevi@easih.ac.uk
#               
#
#   Description:
#   ============
#   script for picking recent file in the input directory - creating fastq files, if need be   - creating md5 files - mounting HPC file system, if not mounted - copying relevant files   from shared disc to HPC - checking md5sum.
#
#*************************************************************************
#
#   Usage:
#   ======
#   $0 -h<help> -q<just create fastq files - nothing else> [input directory(ideal:/data/<project name>/raw), default=current directory]
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  05.05.11 Original.
#   V1.1  06.06.11 Option q added for just creating fastq files, but nothing else.
#   V1.2  07.06.11 Merged fq files are processed instead of fq files.
#*************************************************************************

use strict;
use Getopt::Std;

my %opts;
getopts('hq', \%opts);
#getopts('hpq', \%opts);
#my $copy_files        = $opts{c} || 0;
#my $platform          = $opts{p} || 0;
my $create_fastq_file = $opts{q} || 0;

my $indir = shift;

if($opts{h})
{
  $0 =~ s/.*\///;
  print STDERR "\n\nDescription: Uploads relevant data to HPC\n";
  print STDERR "\nDetails: script for picking recent file in the input directory - creating fastq files, if need be - creating md5 files - mounting HPC file system, if not mounted - copying relevant files from shared disc to HPC - checking md5sum.\n";
  #print STDERR "Usage: $0 -c<copy> -p<platform> -q<create fastq files> -h<help> [input directory(ideally: e.g., /data/<project name>/raw), default = current directory]\n\n";
  print STDERR "\nUsage: $0 -h<help> -q<just create fastq files - nothing else> [input directory(ideal: /data/<project name>/raw), default = current directory]\n\n\n";
  exit;
}

my $pwd;

if($indir)
{
    chomp $indir;
    $indir =~ s/\/$//;
    $pwd = "$indir";
}

chomp($pwd = `pwd`), unless($pwd);
print "\n Input directory (default = current directory): $pwd \n";

my($project) = $1, if($pwd =~ /\/data\/(.*)\/raw/);
print "\n Project name: $project \n", if($project);

my(@files) = PickRecentFile($pwd);

my(@fq_files,@cs_files);

my $flag = "";

foreach(@files)
{
    my $merge = 0;

    if((/(.*\/[A-Z]\d{7}_\d+)(m?).*fq.gz/) || (/(.*\/[A-Z]\d{7}).*fq.gz/))
    {
	$flag = $1;
	$merge = 1, if($2);
	
	if($merge)
	{
	    my $orgline = $_;
	    my $tmpline = $orgline;
	    $tmpline =~ s/m\./\./;
	    print "Not processing fq file: $tmpline, instead processing merged fq file: $orgline \n";	    
	    pop @fq_files;
	    push @fq_files, $orgline;
	}
	else
	{
	    push @fq_files, $_;
	    print "fastq file: $_ \n";
	}
    }
    elsif((/(.*\/[A-Z]\d{7}_\d+).*csfasta.gz/) || (/(.*\/[A-Z]\d{7}).*csfasta.gz/))
    {
	if($1 ne $flag)
	{
	    push @cs_files, $_;
	    print "csfasta file that does not have a corresponding fastq file: $_ \n";
	}	
    }
}

my(@newfq_files) = CreateFastq(@cs_files), if(@cs_files);

exit, if($create_fastq_file);

my(@allfq_files) = sort(@fq_files,@newfq_files);
    
if(@allfq_files)
{
    ComputeMd5sum(@allfq_files);
    UploadToHPC($project,@allfq_files);
}    

    
#################################################################################
sub PickRecentFile
{
    my($pwd) = @_;
    
    chomp(my $test =  `ls -rt ${pwd}/*.gz | tail -1`); 
    print "\n Last updated file in the directory: $test \n";

    chomp(my $recent = `ls -lrt  ${pwd}/*.gz | awk '{print \$6}' | tail -1`); 
    print "\n date of recent file (files which share same date will be pulled out for further processing): $recent \n";

    chomp(my $today = `date '+%F'`);
    
    chomp(my $date_diff = `echo "\$((\$(date -d '$today' "+%s") - \$(date -d '$recent' "+%s")))" | awk '{printf "%.0f", \$1/(24*60*60)}'`);
    
    print "\n\n *** ALERT! Files of interest are older than a day: Date difference is $date_diff day(s) ***\n\n", if(($date_diff>1));

    chomp(my $list = `ls -lrt ${pwd}/*.gz | grep "$recent" | awk '{print \$8}' | sort | xargs  echo`);
    print "\n List of files which share same date as the recent file in the directory: $list \n\n";
    
    my(@files) = split(" ",$list);

    return(@files);
}


#################################################################################
sub CreateFastq
{
    my(@cs_files) = @_;
    my($files,@fq_files);

    
    foreach(@cs_files)    
    {
	chomp;
	
	if(/(.*\/[A-Z]\d{7}.*)_[F|R].*/)
	{
	    my $basename = $1;
	    
	    print "Creating fastq file for: $basename \n";
	    system("~kb468/easih-toolbox/scripts/csfasta2fq.pl -p $basename -o $basename");
	    push @fq_files, "${basename}.1.fq.gz", if(-e "${basename}.1.fq.gz");
	}
    }

    return(@fq_files);
}


#################################################################################
sub UploadToHPC
{
    my($project,@fq_files) = @_;
    my(@md5_files,@hpcfiles);

    foreach(@fq_files)
    {
	my $md5file = "$_.md5";
 	push(@md5_files, $md5file);
    }
    
    my @upload_files = sort(@fq_files,@md5_files);
    	
    my $files = join(" ",@upload_files);

    print "File(s) to upload to HPC: $files \n";


    # Grab the userid 
    chomp(my $userid = `whoami`);
    my $dirhandle = "/home/$userid/hpc/$userid/$project";

    # Check the size of mount point 
    system("mkdir -p /home/$userid/hpc"), unless(-d "/home/$userid/hpc");
    chomp(my $hpc = `ls -l ~$userid | grep hpc | awk '{print \$5}'`); 
                
    print "Mounting HPC filesystem for you ($userid), if it is not currently mounted... \n";

    # Mount size: zero - not mounted; non-zero - mounted
    system("sshfs $userid\@login.hpc.cam.ac.uk:/scratch /home/$userid/hpc"), unless((-d "/home/$userid/hpc") && ($hpc));
    system("mkdir -p $dirhandle"), if(-d "/home/$userid/hpc");
    system("chmod -R g+srwx $dirhandle");
    system("cp $files $dirhandle");

    foreach(@md5_files)
    {
	s/.*([A-Z]\d{7}.*)/${dirhandle}\/$1/;
	open INFILE, "$_" or die $!;
	
	s/(.*gz)(.*)/${1}.hpc${2}/;
	open OUTFILE, ">$_" or die $!;
	push @hpcfiles, $_; 
	
	while(<INFILE>)
	{
	    s/(.*\s+).*([A-Z]\d{7}.*)/${1}${dirhandle}\/${2}/;
	    print OUTFILE "$_";
	}
	
	close INFILE;
	close OUTFILE;
    }

    CheckMd5sum(@hpcfiles);
}


#################################################################################
sub ComputeMd5sum
{

    my(@files) = @_;

    foreach my $file (@files)
    {
	print "Creating md5 file for: $file \n";
	system("md5sum $file > ${file}.md5");
    }
}


#################################################################################
sub CheckMd5sum
{
    my(@files) = @_;

    foreach my $file (@files)
    {
	print "Checking md5sum against HPC file: $file \n";
	chomp(my $not_ok = `md5sum -c $file | grep -v "OK\$"`);
	print "\n\n *** Oops! File integrity check failed - corruption during file transfer?: $not_ok *** \n\n", if($not_ok);

    }
}
