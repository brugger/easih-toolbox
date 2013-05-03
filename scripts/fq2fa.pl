#!/usr/bin/perl


while (<STDIN>){
	print ">$_";
	$dump=<STDIN>;
	print $dump;
	$dump=<STDIN>;
	$dump=<STDIN>;
}
	
