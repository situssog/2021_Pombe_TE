use strict;

##### This script is used to filter out sequences with length < $cutoff.

# Two parameters, example: perl 0.keep_longreads_greaterThan1000bp.pl JB22.longreads.fasta 1000
my $fastafile = $ARGV[0];   # long reads fasta
my $len_cutoff = $ARGV[1];  # read length cutoff


my %all;
my %len;
open(IN,"<$fastafile") || die;
open(OUT,">$fastafile.$len_cutoff.fasta") || die;
my $mark = "nothing";
while (<IN>){
	chomp;
	chop if(/\r$/);
	
	if(/^\>(\S+)/){		
		$mark = $1;
	}else{
		my $tmp = length($_);
		if(defined $all{$mark}){
			$all{$mark} .= $_;
			$len{$mark} += $tmp;
		}else{
			$all{$mark} = $_;
			$len{$mark} = $tmp;
		}
	}
}
close IN;

foreach my $ifasta (keys %all){
	next if($len{$ifasta} < $len_cutoff );
	print OUT ">$ifasta\n$all{$ifasta}\n";
	
}
close OUT;

open(OUT2,">$fastafile.$len_cutoff.lengthreport.txt") || die;
foreach my $ifasta (keys %all){
	next if($len{$ifasta} < $len_cutoff );
	print OUT2 "$ifasta\t$len{$ifasta}\n";
	
}
close OUT2;


