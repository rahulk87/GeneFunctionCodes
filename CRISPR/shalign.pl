#!/usr/bin/perl
use Parallel::ForkManager;
use Getopt::Long;
use warnings;
use strict;
# changes by James: changed location of ccodeFile at line 47; commented 787+788 (fq comment match check); re-wrote 792+793 (regex updated for new format); commented all occurrences of 'unlink' for debugging - this version is now in shalign.pl~ and the unlinks have been restored here for testing without the --dbload option...
sub parseInputFile($$$$);
sub readShrnaLibFile($$);
sub generateBinFromFile();
sub printBinned($$$$);
sub printAmbiguousSummary();
sub doCDistance($$$$$);
sub distance2Bins($$);
sub produceReports($$);
sub printFile($$$);
sub printResults($$$$$$);
sub printAmbiguous($$);
sub printUnmapped($$$$$);
sub printNdist($$$);
sub postprocess($$$);
sub postProcessNdistr($$$);
sub postProcessBinned($$);
sub postProcessUnmapped($$);
sub postProcessDB($$);
sub postProcessResults($$$);
sub printShrnaNgsStats($$$$); # expId, binData, outputfile
sub printShrnaNgsNdist($$$$); # expId, ndistData, outputfile

my $expId = 0;				# Project ID (only used in output bin file)
my $binThres = 0;			# Threshold for bins (if any)
my $offset = 0;			# Left sequence boundary (if any)
my $length = 19;			# Right sequence boundary (if any)
my @fastqFiles = ();			# filenames of input files 
my @shrnaLibFiles = ();			# Name of files containing shrna libraries
my $distanceOutFile = "";		# File to save distances (Hamming or Damerau-Levenshtein)
my $distanceThreshold = 2;		# Threshold for minimum distance (Hamming or Damerau-Levenshtein), set this rather high by default, user should supply lower values
my $srId = 0;				# maximum sr_id value in NGS_SR_SEQ table
my $binned_id_counter = 0;
my %sequences = ();
my @sortedSeqFiles = ();		# sorted fastq files by lane
my %shrnaLibrary = ();
my %ambiguousCalls = ();
my $hammingDistance;
my $damerauLevenshteinDistance;
my $joinReports;
my $outputDirectory;
my $ccodeFile = "/scratch/breakthr/jamesc/scripts/shALIGN/distance";	# Location of the file containing C distance calculation (Hamming and Damerau-Levenshtein

######
my %dbloadStats = (); 	# {lane} ->[sr_total, bin_total, TACA_total, oneN, twoN, sr_perfect_mapped_total,
						# sr_distance1_total, sr_distance2_total, sr_ambiguous, bin_perfect_mapped_total,
						# bin_distance1_total, bin_distance2_total, bin_ambiguous

######

GetOptions (	"project=i"			=> \$expId,
				"h"					=> \$hammingDistance,
				"dl"				=> \$damerauLevenshteinDistance,
				"threshold=i" 		=> \$binThres,
				"fastqFiles=s"		=>	\@fastqFiles,
				"offset=i" 		=> \$offset,
				"length=i"		=> \$length,
				"libraries=s"		=> \@shrnaLibFiles,
				"distance=i"		=> \$distanceThreshold,
				"dbload"		=>\$joinReports,
				"outDir=s"		=> \$outputDirectory,
);

@fastqFiles ? @fastqFiles = split(/,/,join(',', @fastqFiles)) : die "Please provide at least one fastq file.\n";
@shrnaLibFiles = split(/,/,join(',', @shrnaLibFiles));

### Sanity checks
die "Number of library files is different to number of fastq files.\n" if (@shrnaLibFiles > 0 && @shrnaLibFiles != @fastqFiles);
die "please specify length\n" if (! $length);
die "please specify project\n" if ($joinReports && ! $expId);
if ($outputDirectory) {
	if(! $outputDirectory =~ /\/$/) {
		$outputDirectory .= "/";
	}
	open(DIR, $outputDirectory) or die "output directory is not accessible\n";
	close(DIR);
}
###

### TO DO: control thread dispatcher
my $pm = new Parallel::ForkManager(8); # Run in parallel with 8 threads
print "Number of fastq files:\t" . scalar @fastqFiles . "\n";

for (my $i =0; $i < @fastqFiles; $i++) {
	my $fastqFile = $fastqFiles[$i];
	#print STDERR "Processing file $fastqFile\n";
	my $shrnaLibFile;
	@shrnaLibFiles ? $shrnaLibFile = $shrnaLibFiles[$i] : undef $shrnaLibFile;	
	
	# initialise thread
	my $pid = $pm->start($fastqFile) and next;
	my ($binsRef, $npositions) = parseInputFile($fastqFile,$offset, $length, $shrnaLibFile);
	# binsRef structure:
	# {lane}->{trimmed seq}->[hits,TACA,binNs,[distance],[Library hit]]
	

	if ($shrnaLibFile) {
		my $outputFile = $fastqFile . ".distance.out"; 	# this files should not be kept
		my $inputFile = $fastqFile . ".distance.in";	# this file should not be kept
		my $resultsFile = $fastqFile . ".results";
		my $resultsSummaryFile = $resultsFile . ".summary";
		my $ambiguousFile = $fastqFile . ".ambiguous";
		my $shrnaLibrary = readShrnaLibFile($shrnaLibFile, $binsRef);
		printFile($inputFile, $shrnaLibrary, $binsRef); #generate input file for distance computation
		doCDistance($ccodeFile, $inputFile, $shrnaLibFile, $outputFile, $distanceThreshold); # compute distance
		distance2Bins($binsRef, $outputFile);
		my ($resultsData, $ambiguousData, $unmappedData) = produceReports($binsRef, $shrnaLibrary);
		printResults($expId, $shrnaLibrary, $resultsData, $resultsFile, $resultsSummaryFile, $distanceThreshold);
		printAmbiguous($ambiguousData, $ambiguousFile);

		if($joinReports) {
			my $shrnaNgsStatsFh;
			my $shrnaNgsStatsFile = $fastqFile . ".shrnaNgsStats";
			open ($shrnaNgsStatsFh, ">", $shrnaNgsStatsFile) or die "$0: Could not open $shrnaNgsStatsFile $!";
			printShrnaNgsStats($expId, $binsRef, $resultsData,  $shrnaNgsStatsFh);
			close($shrnaNgsStatsFh);
			my $shrnaNgsNdistFh;
			my $shrnaNgsNdistFile = $fastqFile . ".shrnaNgsNdist";
			open ($shrnaNgsNdistFh, ">", $shrnaNgsNdistFile) or die "$0: Could not open $shrnaNgsStatsFile $!";
#			printShrnaNgsNdist($expId, $npositions, $shrnaNgsNdistFh);
			printShrnaNgsNdist($expId, $npositions, $shrnaNgsNdistFh, $length);
			close($shrnaNgsNdistFh);
		}
		
		######## print unmapped & contaminant file
		my $unmappedFh;
		my $unmappedSummaryFh;
		my $contamFh;
		my $unmappedFile = $fastqFile . ".unmapped";
		my $unmappedSummaryFile = $unmappedFile . ".summary";
		my $contaminantsFile = $fastqFile . ".contaminants";
		open($unmappedFh, ">", $unmappedFile) or die "$0: Could not open $unmappedFile $!";
		open($unmappedSummaryFh, ">", $unmappedSummaryFile) or die "$0: Could not open $unmappedSummaryFile $!";
		open($contamFh, ">", $contaminantsFile) or die "$0: Could not open $contaminantsFile $!";
		printUnmapped($expId, $binsRef, $unmappedFh, $unmappedSummaryFh, $contamFh);
		close($unmappedFh);
		close($unmappedSummaryFh);
		close($contamFh);
	}
	### print bin file and summary per lane
	my $binFh;
	my $binSummaryFh;
	my $binHitsProfileFh;
	my $binnedFile = $fastqFile . ".binned";
	my $binFileSummary = $binnedFile . ".summary";
	my $binFileProfile = $fastqFile . ".binprofile";
	open($binFh, ">", $binnedFile) or die "$0: Could not open $binnedFile $!";
	open($binSummaryFh, ">", $binFileSummary) or die "$0: Could not open $binFileSummary $!";
	open($binHitsProfileFh, ">", $binFileProfile) or die "$0: Could not open $binFileProfile $!";
	
	printBinned($binsRef, $binFh, $binSummaryFh, $binHitsProfileFh);
	
	close($binFh);
	close($binSummaryFh);
	close($binHitsProfileFh);
	
	############### print N distribution summary per lane
	my $ndistrFh;
	my $ndistFile = $fastqFile . ".ndistribution";
	open($ndistrFh, ">", $ndistFile) or die "$0: Could not open $ndistFile $!";
	printNdist($npositions, $ndistrFh, $length);
	close($ndistrFh);
	
	#################
	$pm->finish
}
$pm->wait_all_children;
#postprocess(\@fastqFiles, "ndistribution", $ndistrFileTitle);
#postprocess(\@fastqFiles, "binned", $binFileTitle);
if (@shrnaLibFiles) { 
	postProcessResults(\@fastqFiles, "results.summary", $distanceThreshold);
	postProcessUnmapped(\@fastqFiles, "unmapped.summary");
}
postProcessNdistr(\@fastqFiles, "ndistribution", $length);
postProcessBinned(\@fastqFiles, "binned.summary");


if ($joinReports) {
	postProcessDB(\@fastqFiles, "binned");
	postProcessDB(\@fastqFiles, "binprofile");
}

if ($joinReports && @shrnaLibFiles) {
	postProcessDB(\@fastqFiles, "unmapped");
	postProcessDB(\@fastqFiles, "ambiguous");
	postProcessDB(\@fastqFiles, "contaminants");
	postProcessDB(\@fastqFiles, "results");
	postProcessDB(\@fastqFiles, "shrnaNgsStats");
	postProcessDB(\@fastqFiles, "shrnaNgsNdist");

}	

print "Finished:\t" . localtime() . "\n";

#################
# printShrnaNgsNdist($$$)
################
# Print of number of Ns at each position over specified length of
# short reads (first and last cycles)
sub printShrnaNgsNdist($$$$) {
	print STDERR "Writing shrnaNgsNdist file\n";
	my ($expId, $npositions, $ndistrFh, $length) = @_;
	foreach my $lane (sort {$a <=> $b} keys %$npositions) {
		#foreach my $position (sort {$a <=> $b} keys %{$npositions->{$lane}}) {
		for (my $position = 1; $position <= $length; $position++) {
			#print STDOUT $npositions . "\t" . $lane . "\t" . $position . "\n";
			$npositions->{$lane}->{$position} ? (print $ndistrFh "$expId\t$lane\t$position\t" . $npositions->{$lane}->{$position} . "\n") : (print $ndistrFh "$expId\t$lane\t$position\t0\n");
		}
	}
}


##################################################
# sub printShrnaNgsStats($$$$)
##################################################
sub printShrnaNgsStats($$$$) {
	print STDERR "Writing shrnaNgsStatsTable data\n";
	my ($expId, $binsRef, $resultsData, $ngsStatsFh) = @_;
	my %binThres = ();

	foreach my $lane (sort keys %$binsRef) {
		my $TACAsum = 0; my $n1sum = 0; my $n2sum = 0; my $nsum = 0; my $srsum = 0;
		my $binsum = 0; my $unmappedsum = 0; my $unmappedbinsum = 0; my $perfectsum = 0;
		my $perfectbinsum = 0; my $ambiguoussum = 0; my $ambiguousbinsum = 0;
		my $dist1sum = 0; my $dist1binsum = 0; my $dist2sum = 0; my $dist2binsum = 0;
		
		# compute the inexact match sums from resutlsData
		foreach my $shrna (keys %{$resultsData->{$lane}}) {
			$resultsData->{$lane}->{$shrna}->[1] ? (($dist1sum += $resultsData->{$lane}->{$shrna}->[1] ) && $dist1binsum++) : 0;
			$resultsData->{$lane}->{$shrna}->[2] ? (($dist2sum += $resultsData->{$lane}->{$shrna}->[2] ) && $dist2binsum++) : 0;
		}
		foreach my $shortRead (keys %{$binsRef->{$lane}}) {
			$binsum++;
			$srsum += $binsRef->{$lane}->{$shortRead}->[0];
			if (exists $binsRef->{$lane}->{$shortRead}->[5]) {
				if ($binsRef->{$lane}->{$shortRead}->[5] eq "mapped") {
					$perfectbinsum++;
					$perfectsum += $binsRef->{$lane}->{$shortRead}->[0];
				}
				elsif ($binsRef->{$lane}->{$shortRead}->[5] eq "unmapped") {
					$unmappedbinsum++;
					$unmappedsum += $binsRef->{$lane}->{$shortRead}->[0];
				}
				elsif ($binsRef->{$lane}->{$shortRead}->[5] eq "ambiguous") {
					$ambiguousbinsum++;
					$ambiguoussum += $binsRef->{$lane}->{$shortRead}->[0];
				}
			}
			else {		
				$unmappedbinsum++;
				$unmappedsum += $binsRef->{$lane}->{$shortRead}->[0];
			}

			$binsRef->{$lane}->{$shortRead}->[1] ? ($TACAsum += $binsRef->{$lane}->{$shortRead}->[1] ) : 0;# TACAS
			$binsRef->{$lane}->{$shortRead}->[6] ? ($n1sum += $binsRef->{$lane}->{$shortRead}->[6] ) : 0; # n1s
			$binsRef->{$lane}->{$shortRead}->[7] ? ($n2sum += $binsRef->{$lane}->{$shortRead}->[7] ) : 0; # n2s
		}
		print $ngsStatsFh "$expId\t$lane\t$srsum\t$binsum\t$TACAsum\t$n1sum\t$n2sum\t$perfectsum\t$dist1sum\t$dist2sum\t$ambiguoussum\t" . 
								"$unmappedsum\t$perfectbinsum\t$dist1binsum\t$dist2binsum\t$ambiguousbinsum\t$unmappedbinsum\n";			
	}
} 

sub postProcessDB($$) {
	my ($fastqFiles, $identifier) = @_;
	my $outputFile = "report" . ".$identifier";
	local $/;
	open(my $fh, ">", $outputFile) or die "$0: Could not open $outputFile $!";
	foreach my $file (sort @$fastqFiles) {
		my $inputFile = $file . ".$identifier";
		open(IN, $inputFile) or die "$0: Could not open $inputFile $!";
		print $fh <IN>;
		unlink($inputFile);
	}
	close($fh);
}
sub postProcessUnmapped($$) {
my ($fastqFiles, $identifier) = @_;
	my $outputFile = "report" . ".$identifier";
	my @results = ();
	my @sums = ();
	my $outputString = "Lane\tSRs\tbins\n";
	open(my $fh, ">", $outputFile) or die "$0: Could not open $outputFile $!";
	foreach my $file (sort @$fastqFiles) {
		my $inputFile = $file . ".$identifier";
		open(IN, $inputFile) or die "$0: Could not open $inputFile $!";
		while(defined(my $line = <IN>)) {
			chomp($line);
			push(@results, $line);
		}
		close(IN);
		unlink($inputFile) or die "Can't unlink $inputFile: $!";
	}
	foreach my $result (sort @results) {
		$outputString .=  $result . "\n";
		my @fields = split("\t", $result);
		for (my $i = 1; $i < scalar (@fields); $i++) {
			$sums[$i] += $fields[$i];
		}
	}
	$outputString .= "Total";
	for (my $j = 1; $j < scalar (@sums); $j++) {
		$outputString .= "\t" . $sums[$j];
	}
	$outputString .= "\n";
	print $fh $outputString;
	close($fh);
}


sub postProcessBinned($$) {
	my ($fastqFiles, $identifier) = @_;
	my $outputFile = "report" . ".$identifier";
	my @results = ();
	my @sums = ();
	my $outputString = "Lane\tSRs\tbins\tTACA ends\tN calls\tn1 calls\tn2 calls\n";
	open(my $fh, ">", $outputFile) or die "$0: Could not open $outputFile $!";
	foreach my $file (sort @$fastqFiles) {
		my $inputFile = $file . ".$identifier";
		open(IN, $inputFile) or die "$0: Could not open $inputFile $!";
		while(defined(my $line = <IN>)) {
			chomp($line);
			push(@results, $line);
		}
		close(IN);
		unlink($inputFile) or die "Can't unlink $inputFile: $!";
	}
	foreach my $result (sort @results) {
		$outputString .=  $result . "\n";
		my @fields = split("\t", $result);
		for (my $i = 1; $i < scalar (@fields); $i++) {
			$sums[$i] += $fields[$i];
		}
	}
	$outputString .= "Total";
	for (my $j = 1; $j < scalar (@sums); $j++) {
		$outputString .= "\t" . $sums[$j];
	}
	$outputString .= "\n";
	print $fh $outputString;
	close($fh);
}


sub postProcessNdistr($$$) {
	my ($fastqFiles, $identifier, $length) = @_;
	my $outputFile = "report" . ".$identifier";
	my @results = ();
	my @sums = ();
	my $outputString = "Lane";
	for (my $column = 1; $column <= $length; $column++) {
		$outputString .= "\tP$column";
	}
	$outputString .= "\tTotal\n";
	open(my $fh, ">", $outputFile) or die "$0: Could not open $outputFile $!";
	foreach my $file (sort @$fastqFiles) {
		my $inputFile = $file . ".$identifier";
		open(IN, $inputFile) or die "$0: Could not open $inputFile $!";
		while(defined(my $line = <IN>)) {
			chomp($line);
			push(@results, $line);
		}
		close(IN);
		unlink($inputFile) or die "Can't unlink $inputFile: $!";
	}
	foreach my $result (sort @results) {
		my @fields = split("\t", $result);
		$outputString .= $result . "\n";
		for (my $i = 1; $i < scalar (@fields); $i++) {
			$sums[$i] += $fields[$i];
		}
	}
	$outputString .= "Total";
	for (my $j = 1; $j < scalar (@sums); $j++) {
		$outputString .= "\t" . $sums[$j];
	}
	$outputString .= "\n";
	print $fh $outputString;
	close($fh);

}

#################
# postProcessResults($$$)
################
sub postProcessResults($$$) {
	my ($fastqFiles, $identifier, $distanceThreshold) = @_;
	my $outputFile = "report" . ".$identifier";
	open(my $fh, ">", $outputFile) or die "$0: Could not open $outputFile $!";
	my @results = ();
	my @sums = ();
	# print header
	my $outputString = "Lane";
	for (my $i = 0; $i <= $distanceThreshold; $i++) {
		$outputString .= "\tD[$i]";
	}
	$outputString .= "\tTotal\n";
	foreach my $file (sort @$fastqFiles) {
		my $inputFile = $file . ".$identifier";
		open(IN, $inputFile) or die "$0: Could not open $inputFile $!";
		while(defined(my $line = <IN>)) {
			chomp($line);
			push(@results, $line);
		}
		close(IN);
		unlink($inputFile) or die "Can't unlink $inputFile: $!";
	}
	foreach my $result (sort @results) {
		my @fields = split("\t", $result);
		$outputString .= $result . "\n";
	}
	print $fh $outputString;
	close($fh);
}
				
#################
# printNdist($$)
################
# Print of number of Ns at each position over specified length of
# short reads (first and last cycles)
sub printNdist($$$) {
	print STDERR "Writing N distribution file\n";
	my ($npositions, $ndistrFh, $length) = @_;
	foreach my $lane (sort {$a <=> $b} keys %$npositions) {
		my $laneTotal = 0;
		#my ($sequenceLength) = scalar (keys %{$npositions->{$lane}});
		#$title .=  "Lane";
		#for (my $i=1; $i <= $sequenceLength; $i++) {
	#		$title .= "\t[$i]";
	#	}
	#	$title .="\n";
		print $ndistrFh "$lane";
		#foreach my $position (sort {$a <=> $b} keys %{$npositions->{$lane}}) {
		for (my $position = 1; $position <= $length; $position++) {
			#print STDOUT $npositions . "\t" . $lane . "\t" . $position . "\n";
			$npositions->{$lane}->{$position} ? ((print $ndistrFh "\t" . $npositions->{$lane}->{$position}) && ($laneTotal +=$npositions->{$lane}->{$position})) : print $ndistrFh "\t0";
		}
		print $ndistrFh "\t$laneTotal\n";
	}
}

##################################################
# sub printBinned($$$$)
##################################################
sub printBinned($$$$) {
	print STDERR "Writing binned data\n";
	my ($binsRef, $binFh, $binSummaryFh, $binHitsProfileFh) = @_;
	my %binThres = ();
	if (! $joinReports) { print $binFh "Lane\tSR sequence\tbin hits\tclass\tTACA endings\tN calls\tn1 calls\tn2 calls\n";}
	foreach my $lane (sort keys %$binsRef) {
		my $TACAsum = 0;
		my $n1sum = 0;
		my $n2sum = 0;
		my $nsum = 0;
		my $srsum = 0;
		my $binsum = 0;
		foreach my $shortRead (keys %{$binsRef->{$lane}}) {
			print $binFh "$lane\t$shortRead\t" . $binsRef->{$lane}->{$shortRead}->[0] . "\t"; # hits;
			$binsRef->{$lane}->{$shortRead}->[5] ? (print $binFh $binsRef->{$lane}->{$shortRead}->[5] . "\t") : 
					print $binFh "unmapped\t"; # match status
			$binThres{$lane}->{ $binsRef->{$lane}->{$shortRead}->[0] }++;
			$binsRef->{$lane}->{$shortRead}->[1] ? ((print $binFh $binsRef->{$lane}->{$shortRead}->[1] . "\t") && ($TACAsum += $binsRef->{$lane}->{$shortRead}->[1] )) : # TACAS
					print $binFh "0\t";
			$binsRef->{$lane}->{$shortRead}->[2] ? ((print $binFh $binsRef->{$lane}->{$shortRead}->[2] . "\t") && ($nsum += $binsRef->{$lane}->{$shortRead}->[2])) : # Ns in bins
					print $binFh "0\t";
			$binsRef->{$lane}->{$shortRead}->[6] ? ((print $binFh $binsRef->{$lane}->{$shortRead}->[6] . "\t") && ($n1sum += $binsRef->{$lane}->{$shortRead}->[6] )) : # n1s
					print $binFh "0\t";	
			$binsRef->{$lane}->{$shortRead}->[7] ? ((print $binFh $binsRef->{$lane}->{$shortRead}->[7] . "\n") && ($n2sum += $binsRef->{$lane}->{$shortRead}->[7] )) : # n2s
					print $binFh "0\n";						
			$binsum++;
			$srsum += $binsRef->{$lane}->{$shortRead}->[0];
		}
		print $binSummaryFh "$lane\t$srsum\t$binsum\t$TACAsum\t$nsum\t$n1sum\t$n2sum\n";
		if (! $joinReports) { print $binHitsProfileFh "Lane\tbin size\tbin frequency\n";}
		foreach my $lane (sort keys %binThres) {
			foreach my $binsize (sort {$a <=> $b} keys %{$binThres{$lane}}) {
				if ($joinReports && $expId > 0) {
					print $binHitsProfileFh "$expId\t";
				}
				print $binHitsProfileFh "$lane\t$binsize\t" . $binThres{$lane}->{$binsize} . "\n";
			}
		}			
	}
}
##################################################
# sub printUnmapped($$$$)
##################################################
sub printUnmapped($$$$$) {
	print STDERR "Writing unmapped hits to file\n";
	my ($expId, $binsRef, $unmappedFh, $unmappedSummaryFh, $contaminantsFh) = @_;
	if (! $joinReports) { 
		print $unmappedFh "Lane\tSR sequence\tbin hits\tstatus\n";
		print $contaminantsFh "Lane\tSR sequence\tbin hits\n";
	}
	my $count50 = 50;
	foreach my $lane (sort keys %$binsRef) {
		my $srSum = 0;
		my $binSum = 0;
		foreach my $binSequence (sort {$binsRef->{$lane}->{$b}->[0] <=> $binsRef->{$lane}->{$a}->[0]} keys %{$binsRef->{$lane}}) {
			if ($binsRef->{$lane}->{$binSequence}->[5] eq "unmapped") {
				$count50--;
				$binSum++;
				$srSum += $binsRef->{$lane}->{$binSequence}->[0];
				($count50 >= 0) ? (($expId > 0 ? print $contaminantsFh "$expId\t" : 0) && (print $contaminantsFh "$lane\t$binSequence\t" . $binsRef->{$lane}->{$binSequence}->[0] . "\n")) : 0;
				print $unmappedFh "$lane\t$binSequence\t" . $binsRef->{$lane}->{$binSequence}->[0] . "\t" .
					$binsRef->{$lane}->{$binSequence}->[5] . "\n";
			}
		}
		print $unmappedSummaryFh "$lane\t$srSum\t$binSum\n";
	}
}

##################################################
# sub printAmbiguous($$)
##################################################
sub printAmbiguous($$) {
	print STDERR "Writing ambiguous distance hits to file\n";
	my ($ambiguousHits, $outputFile) = @_;
	open(AMB, ">", $outputFile) or die "Could not create $outputFile ambiguous hits file\n";
	if (! $joinReports) {print AMB "Lane\tbinned sequence\tdistance\tshrna\n";}
	foreach my $lane ( sort keys %$ambiguousHits ) {
		foreach my $binSequence (keys %{$ambiguousHits->{$lane}}) {
			foreach my $ambiguousHit (@{$ambiguousHits->{$lane}->{$binSequence}}) {
				print AMB "$lane\t$binSequence\t$ambiguousHit\n";
			}
		}
	}
	close (AMB) or die "Could not close file $outputFile\n";
}
##################################################
# sub printResults($$$$$)
# prints a summary of hits for each shrna
# at each distance as well as the total number
# of hits and the ratio of perfect matches
# to total matches
##################################################
sub printResults($$$$$$) {
	print STDERR "Writing results\n";
	my ($expId, $shrnaLibrary, $resultsData, $outputFile, $resultsSummaryFile, $distanceThreshold) = @_;
	my %summaryStats=();
	open (RES, ">", $outputFile) or die "Could not create $outputFile results file\n";
	open (RESSUM, ">", $resultsSummaryFile) or die "Could not create $resultsSummaryFile $!";
	# print header
	if (! $joinReports) {
		print RES "Lane\tshrna id\tshrna sequence";
		for (my $index = 0; $index <= $distanceThreshold; $index++) {
			print RES "\tdistance=" . $index;
		}
		print RES "\ttotal hits\tmatched ratio\n";
	}
	foreach my $lane (sort keys %$resultsData) {
		foreach my $shrna (keys %$shrnaLibrary) {
			my $shrnaHits = 0;
			# print lane, shrna id and shrna
			my $string;
			# print experiment id if we are in dbload mode
			$expId > 0 ? ($string .= $expId . "\t") : 0;
			$string .= "$lane\t" . ($shrnaLibrary->{$shrna}) . "\t$shrna" ;
			
			# print hits for  distances equal or greater than zero and increment total
			# hits counter accordingly. Otherwise print zero hits
			for (my $i = 0; $i <= $distanceThreshold; $i++) {
				 ($resultsData->{$lane}->{$shrna}->[$i]) ?
					(($string .= "\t" . $resultsData->{$lane}->{$shrna}->[$i]) 
						&& ($shrnaHits += $resultsData->{$lane}->{$shrna}->[$i])
						&& ($summaryStats{$lane}->[$i]->[0]++)
						&& ($summaryStats{$lane}->[$i]->[1] += $resultsData->{$lane}->{$shrna}->[$i])) :
							($string .= "\t0");		
			}
			# print total hits
			$string .= "\t$shrnaHits";
			# print perfect matches to total hits ratio unless we are in dbload mode
			if ($expId == 0) {
				($resultsData->{$lane}->{$shrna}->[0] && $shrnaHits > 0) ? 
					($string .=  "\t" . 
						(sprintf "%4.2f", ($resultsData->{$lane}->{$shrna}->[0] / $shrnaHits))) :
							($string .= "\t0.00");
			}
			$string .= "\n";
			print RES $string; # if (0 < $shrnaHits); ## CAN BE CONTROLLED TO PRINT ALL
		}
	}		
	close(RES) or die "Could not close file $outputFile\n";
	#### print results summary file for the lane
	my $outputString = "";
	foreach my $lane (sort keys %summaryStats) {
		$outputString .= $lane;
		my $srSum = 0;
		my $binSum = 0;
		foreach my $distance (@{$summaryStats{$lane}}) {
			$srSum += $distance->[1];
			$binSum += $distance->[0];
		}
		foreach my $distance (@{$summaryStats{$lane}}) {
			$outputString .= "\t" . (sprintf "%4.2f", 100 * $distance->[1] / $srSum) . "% (" . (sprintf "%4.2f", 100 * $distance->[0] / $binSum) . "%)";
		}
		$outputString .=  "\t$srSum ($binSum)\n";
		print RESSUM $outputString;
	}
	close(RESSUM) or die "Could not close file $resultsSummaryFile\n";
}


###########
# Produce reports code needs rectifying!!!
#########
sub produceReports($$) {
	print STDERR "Producing reports\n";
	my ($binsRef, $shrnaLibrary) = @_;
	my %resultsData = ();
	my %ambiguousData = ();
	my %unmappedData = ();
	foreach my $lane (keys %$binsRef) {
		foreach my $binSequence (keys %{$binsRef->{$lane}}) {
			# add direct library hits to results
			if ( $shrnaLibrary->{$binSequence} ) {
				$resultsData{$lane}->{$binSequence}->[0] = $binsRef->{$lane}->{$binSequence}->[0];
			}
			# investigate indirect hits:
			if ( $binsRef->{$lane}->{$binSequence}->[3] ) {
				my @distances = @{$binsRef->{$lane}->{$binSequence}->[3]};
				if (not @distances) { # no distance mapping
					$unmappedData{$lane}->{$binSequence} = $binsRef->{$lane}->{$binSequence};
					$binsRef->{$lane}->{$binSequence}->[5] = "unmapped";
				}
				else {
					my %hash = ();
					foreach my $distance (@distances) {
						$hash{$distance}++;
					}
					my @sorted = (sort keys %hash);
					if ($hash{$sorted[0]} > 1) {
						for (my $i = 0; $i < @distances; $i++) {
							if ($distances[$i] == $sorted[0]) { 
								push (@{$ambiguousData{$lane}->{$binSequence}},  ($distances[$i] . "\t" . 
									$binsRef->{$lane}->{$binSequence}->[4]->[$i]));
							}
						}
						$binsRef->{$lane}->{$binSequence}->[5] = "ambiguous";	
					}
					else {
						# print to results
						my $libSeq = $binsRef->{$lane}->{$binSequence}->[4]->[0];
						my $dist = $binsRef->{$lane}->{$binSequence}->[3]->[0];
						my $mismatchHits = $binsRef->{$lane}->{$binSequence}->[0];
						$resultsData{$lane}->{$libSeq}->[$dist] += $binsRef->{$lane}->{$binSequence}->[0];
						$binsRef->{$lane}->{$binSequence}->[5] = "distance";
					}
				}
			}
		}
	}
	return (\%resultsData, \%ambiguousData, \%unmappedData);
}

##########
# distance2Bins
##########
sub distance2Bins($$) {
	print STDERR "Doing distance2Bins\n";
	my ($binsRef, $outputFile) = @_;
	open(DIST, $outputFile) or die "Could not read distance file $outputFile\n";
	while(defined(my $line = <DIST>)) {
		#print STDERR $line;
		chomp($line);
		# values[0] -> lane
		# values[1] -> binned sequence
		# values[2] -> bin hits
		# values[3] -> library sequence
		# values[4] -> library seq identifier
		# values[5] -> distance
		my @values = split(/\t/, $line);
		push (@{$binsRef->{$values[0]}->{$values[1]}->[3]}, $values[5]); # distance
		push (@{$binsRef->{$values[0]}->{$values[1]}->[4]}, $values[3]); # shrna
	}
	close(DIST) or die "Could not close file $outputFile\n";
	unlink($outputFile);
}

##########
# prints the levenshtein in file
###########
sub printFile($$$) {
	print STDERR "Writting file for distance compute\n";
	my ($file, $shrnaLibrary, $binsRef) = @_;
	my $hitscounter = 0;
	open(LEV, ">", $file)	
				or die "Could not open $file for writing.\n";
				
	foreach my $lane (keys %$binsRef) {
		foreach my $sr (keys %{$binsRef->{$lane}}) {		
			if(! exists $shrnaLibrary->{$sr}) {
				$binsRef->{$lane}->{$sr}->[5] = "unmapped";
				print LEV "$lane\t$sr\t" . $binsRef->{$lane}->{$sr}->[0] . "\n";
			}
			else {
				$binsRef->{$lane}->{$sr}->[5] = "mapped";
				$hitscounter++;
			}
		}
	}	
	close(LEV) or die "Could not close $file\n";
	print "Matches:\t$hitscounter\n";
}		

###########################
# sub readShrnaLibFile($$)
##########################
sub readShrnaLibFile($$) {
	print STDERR "Reading library file\n";
	my ($shrnaLibFile, $binsRef) = @_;
	my %shrnaLibrary=();
	my @lanes = (keys %$binsRef); 
	open (SHRNA, $shrnaLibFile) or die "Could not read $shrnaLibFile.\n";
		
	while(defined(my $line = <SHRNA>)) {
		chomp($line);
		my @values = split(/\t/, $line);
		$shrnaLibrary{$values[1]} = $values[0];
	}
	close(SHRNA) or die "Could not close file $shrnaLibFile\n";
	return \%shrnaLibrary;
}
###########################
#########################
#               	#
# doCDistance($$$$$)	#
#               	#
#########################

sub doCDistance($$$$$) {
	my ($ccodeFile) = shift @_;
	my ($binnedFile) = shift @_;
        my ($libFile) = shift @_;
        my ($outputFile) = shift @_;
        my ($threshold) = shift @_;

	my @args = ($ccodeFile, $binnedFile, $libFile, $outputFile, $threshold);
	print "$ccodeFile\t$binnedFile\t$libFile\t$outputFile\t$threshold\n";
	my $distanceResult = system(@args);
	print "Finished $binnedFile with exit code\t" . $distanceResult . "\n";
	unlink($binnedFile);
}	

	
#########################
#                  	#
# parseInputFile($)	#
#                  	#
#########################
#
# The script sees that the fastq header 
# contains information about the lane of the run
# If the header is not formated as lane:tile:x:y 
# it complains and exits
#
sub parseInputFile($$$$) {
	print STDERR "Parsing fastq file\n";
	my ($fastqFile, $offset, $length, $shrnaLibFile) = @_;
	my %bins = ();
	my %mismatchBias = ();
	open(IN, $fastqFile) or die "Could not access file $fastqFile\n";
	while(defined(my $line = <IN>)) {
		# Read one fastq record at a time (four lines).
		# Each fastq record comprises four lines
		# Line 1 is the sequence header starting with '@'
		# Line 2 is the short read sequence
		# Line 3 is the phred quality header (similar to line 1)
		# Line 4 is the phred quality score
		if ($line =~ /^@/)	{
			$srId++;
			chomp($line);
			my $seqHeader = substr($line, 1);	# sequence header offset one
			$line = <IN>; #	short read sequence
			chomp($line);
			my $sequence = $line;
			$line = <IN>;
			chomp($line);
			my $phredHeader = substr($line, 1); # phred header offset one
			$line = <IN>;	# phred score (unused)
			chomp($line);
			# perform more sanity checks
#			die "Invalid fastq record in $fastqFile:\tsequence header ( $seqHeader ) different to quality header ( $phredHeader )\n" 
#				if ($seqHeader ne $phredHeader);
			die "Invalid fastq file $fastqFile:\tsequence has different length to quality score for record $seqHeader\n" 
				if (length($sequence) != length($line) );
			# Extract lane from sequence header.  If not available use sample name or exit if not available.
			# modified 12-04-12 by JamesC to look for new, then old style fastq headers
                        my $lane = undef;
                        if($seqHeader =~ /[^:]+:[^:]+:[^:]+:(\d+):[^:]+:[^:]+:[^ ]+ \d+:[NYny]:[^:]+:/){
                          $lane = $1;
			}
                        elsif($seqHeader =~ /[^:]+:(\d):\d+:-?\d+:-?\d+/){
                          $lane = $1;
                          warn("Using older fastq header format to find lane information.\n");
                        }
                        else{
                          $lane = 'proton';
			  #warn("Using proton in place of lane number\n");
			  # die "Inappropriate fastq Header ( $seqHeader ) in file $fastqFile:  Could not deduce lane\n";
                        }
			# If no length for short reads have been specified
			# set it up so that the whole sequence string is used:
			if (0 == $length) { $length = length($sequence); }
			#substring sequence to range of interest (with array offsetting):
			my $seedSequence = substr($sequence, $offset, $length);
			# increment hit counter
			$bins{$lane}->{$seedSequence}->[0]++;

			# register a TACA hit
			if ($sequence =~ /TACA$/) {
				$bins{$lane}->{$seedSequence}->[1]++;
			}
			# store position of Ns for each lane
			# and number of Ns for each bin
			my $nCounter = 0;
			while ($seedSequence =~ /N/g) {
				$nCounter++;
				$bins{$lane}->{$seedSequence}->[2]++;
				$mismatchBias{$lane}->{pos($seedSequence)}++;
			}
			if ($nCounter == 1) {
				$bins{$lane}->{$seedSequence}->[6]++; # 1 N number per bin
			}
			if ($nCounter == 2) {
				$bins{$lane}->{$seedSequence}->[7]++; # 2 N number per bin
			}
			
		}
	}
	close(IN) or die "Could not close file $fastqFile\n";
	print "Processed file $fastqFile\t" . localtime() . "\n";
	return (\%bins, \%mismatchBias);
}
		
