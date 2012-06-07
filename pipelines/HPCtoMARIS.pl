#!/usr/bin/perl -w
use strict;

my @basenames = @ARGV;
my $bedfile = "/data/SureSelect_All_Exon_50mb_with_annotation.hg19.bed";

chomp(my $userid = `whoami`);

foreach my $basename(@basenames)
{
    my(@files,$project);

    $project = $1, if($basename =~ /^([A-Z]\d{2})/);
          
    chomp(my $filelist = `ssh login.hpc.cam.ac.uk "ls scratch/$project/$basename*"`);
 
    @files = split("\n", $filelist);

    foreach my $file(@files)
    {
	next, if($file =~ /\.maris$|\.out$|\.md5$|\.log$|\.gz/);

	system("mkdir -p /data/$project/MARIS");
	system("scp login.hpc.cam.ac.uk:/home/$userid/$file /data/$project/MARIS"),  
    }

    system("cd /data/$project/MARIS; perl /home/svvd2/easih-toolbox/scripts/Variation_report.pl -B $bedfile -b ${basename}.bam -Q $basename"); 
}
  



