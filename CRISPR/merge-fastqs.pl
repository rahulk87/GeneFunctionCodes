#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $fastqs_to_merge;
GetOptions (
  "fastqs_to_merge=s" => \$fastqs_to_merge
  );

my @fastq_files = split(/,/,$fastqs_to_merge);

my $fastq_dir = $fastq_files[0];
$fastq_dir =~ s/[^\/]+$//;

open OUT, "> $fastq_dir/merged.fastq" or die "Unable to write merged.fastq: $!\n";

foreach my $fastq_file (@fastq_files){
	open IN, "< $fastq_file" or die "Can't read input fastq $fastq_file: $!\n";
	while(<IN>){
	  if(/^\@HISEQ/){
		my $header = $_;
		$header =~ s/(\@HISEQ:[^:]+:[^:]+):\d+:/$1:0:/;
			# @HISEQ:69:H2MWKADXY:2: becomes
			# @HISEQ:69:H2MWKADXY:0:
		print OUT $header;
	  }
	  else{
		print OUT $_;
	  }
	}
	close IN;
}

close OUT;