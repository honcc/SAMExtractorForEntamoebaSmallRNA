#!/usr/bin/perl/ -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Getopt::Long;

######################################################################################################################################################
#
#	Description
#		This is a perl script to extract line in a SAM file specifically designed to pickup mismatched small RNA with mismatched polyA at the end.
#
#	Input
#		--bamPath=					a bam file, which is have to be sorted;
#		--hardCodedCriteria=		"yes" or "no". If "yes", the reads will be selected by some extra hard coded criteria as coded; The aim to is give ad hoc flexibility for selection;
#		--misMatchMode=				'post26ntATail' or 'freeLenFreeTail' or 'freeLenATail'; [post26ntATail] post26ntATail pick up reads with most A mismatch after 26nt, freeLenFreeTail pick up reads of any length with certain number of mismatches in the tail;
#		--countMismatchOnly= 		"yes" or "no"; [no] count mismatch only, not print seq; default = no
#		--removeRedundant=			"yes" or "no"; [no]; To remove redundant read or not; redundant read is defined by e.g. "${${${${$cntgPosLenHsh{$curntCntg}}{$readStart}}{$readStrand}}{$readLength}}{"noTail"}++", "noTail" refers to a category without tail;
#		--pileupCounterPath=		path of pileupCounter; it will run pileup counter for the 4 bams; 
#		--GFFPath=					path of Gff for pileupCounter;
#		--fastaPath=				path of fasta for pileupCounter;
#		--outDir=					directory for output, default = a msmatchExtracted folder in the dir of SAM path;
#
#	Usage
#		perl SAMExtractorForEntamoebaSmallRNA_v0.1.pl --misMatchMode=freeLenFreeTail --bamPath=/Volumes/B_MPro2TB/NGS/results/HM1SmallRNA/dsGFP_combined_all_28To35nt/1208_HM1_dsGFP_TPooled_TAPminus/dsGFP_combined_all_28To35nt.1208_HM1_dsGFP_TPooled_TAPminus.sorted.sam --outDir=/Volumes/B_MPro2TB/NGS/results/HM1SmallRNA/dsGFP_combined_all_28To35nt/1208_HM1_dsGFP_TPooled_TAPminus/SAMExtractorForEntamoebaSmallRNA/ --hardCodedCriteria=yes
#
#	Assumption
#
#	History:
#		
#		V0.1
#			-debut, based on SAMExtractor_v0.8
#
#		V0.2
#			-bamPath now takes bam file also;
#			-the output will be piped to samtools directly so all the intermediate sam file outputs will be skipped;
#			-pileupCounterPath option added;
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#
#1----------Read parameters ----------#
printCMDLogOrFinishMessage("CMDLog");

my ($bamPath, $hardCodedCriteria, $misMatchMode, $countMismatchOnly, $removeRedundant, $pileupCounterPath, $GFFPath, $fastaPath, $outDir) = readParameters();

my ($outBamPathHsh_ref, $misMatchCountHsh_ref, $redundantPosHsh_ref) = readAndEditSAMOnTheFly($bamPath, $hardCodedCriteria, $misMatchMode, $countMismatchOnly, $removeRedundant, $outDir);

printMisMatchCount($misMatchCountHsh_ref, $misMatchMode, $outDir);

printRedundantPos($redundantPosHsh_ref, $outDir);

runSamtoolIndexBAMAndPileupCounter($outBamPathHsh_ref, $pileupCounterPath, $GFFPath, $fastaPath);

printCMDLogOrFinishMessage("finishMessage");

exit;

#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	my ($bamPath, $hardCodedCriteria, $misMatchMode, $countMismatchOnly, $removeRedundant, $pileupCounterPath, $GFFPath, $fastaPath, $outDir);
	
	$hardCodedCriteria = "yes";
	$misMatchMode = 'post26ntATail';
	$countMismatchOnly = 'no';
	$removeRedundant = 'no';
	$outDir = 'extractedBAM';
	
	GetOptions (
					"bamPath=s" => \$bamPath,
					"hardCodedCriteria:s" => \$hardCodedCriteria,
					"misMatchMode:s" => \$misMatchMode,
					"countMismatchOnly:s" => \$countMismatchOnly,
					"removeRedundant:s" => \$removeRedundant,
					"pileupCounterPath=s" => \$pileupCounterPath,
					"GFFPath=s" => \$GFFPath,
					"fastaPath=s" => \$fastaPath,
					"outDir:s" => \$outDir,

	) or die "Error in command line arguments\n";
	
	#---check the files
	open (TEST, "$bamPath") || die "Can't open $bamPath\n"; close TEST;
	open (TEST, "$pileupCounterPath") || die "Can't open $pileupCounterPath\n"; close TEST;
	open (TEST, "$GFFPath") || die "Can't open $GFFPath\n"; close TEST;
	open (TEST, "$fastaPath") || die "Can't open $fastaPath\n"; close TEST;

	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	system "mkdir -p -m 777 $outDir/log/pileupCounterRunlog/";
	
	return ($bamPath, $hardCodedCriteria, $misMatchMode, $countMismatchOnly, $removeRedundant, $pileupCounterPath, $GFFPath, $fastaPath, $outDir);
}
########################################################################## readAndEditSAMOnTheFly
sub readAndEditSAMOnTheFly {

	my ($bamPath, $hardCodedCriteria, $misMatchMode, $countMismatchOnly, $removeRedundant, $outDir) = @_;
	
	#---read the SAM
	my ($fileTotalLineNum,  $intervalSize);
	
	if ($bamPath =~ m/\.bam$/) {#----bam file as input
		open (INBAM, " samtools view $bamPath | ");
		($fileTotalLineNum,  $intervalSize) = checkBAMAlignNum($bamPath);
	} else {
		die "bamPath have to ended in either .bam. Quitting\n";
	}
	
	#---define the variables
	my %cntgPosHsh;
	my %redundantPosHsh;

	my @bamPathSplt = split /\//, $bamPath;
	my %outBamPathHsh;
	my %samFileHandleHsh;
	my $sRNACtgryCountHsh_ref = {};

	#---define the start time and counters
	my $intervalStart = time();
	my $lineProc = my $progCount = my $totalTimeSpent = 0;

	print "Start process the alignments.\n";
	printProgressScale("Extracting the BAM file", 50);

	#---Start reading bamPath
	my %misMatchCountHsh;
	
	my ($setCutoffLen, $minPostCutoffLenAMismatchPrptn, $minMisMatchNumPostCutoffLen, $minPostCutoffLenASeqPrptn);

	if ($misMatchMode eq "post26ntATail") {
		$setCutoffLen = "post26ntATail";
		$minPostCutoffLenAMismatchPrptn = 0.3;#---proportion of seq have to be mismatch after the cuttoffLength
		$minPostCutoffLenASeqPrptn = 0.7;#---proportion of seq have to be A after the cuttoffLength
		$minMisMatchNumPostCutoffLen = 2;#----minimum number of mismatches after cuttoffLength
	} elsif ($misMatchMode eq "freeLenFreeTail") {
		$setCutoffLen = "freeLen";
		$minPostCutoffLenAMismatchPrptn = 0;#---proportion of seq have to be mismatch after the cuttoffLength
		$minPostCutoffLenASeqPrptn = 0; #---proportion of seq have to be A after the cuttoffLength
		$minMisMatchNumPostCutoffLen = 3;#----specify the minimum number of mismatches after cutoffLen
	}  elsif ($misMatchMode eq "freeLenATail") {
		$setCutoffLen = "freeLen";
		$minPostCutoffLenAMismatchPrptn = 0.3; #---proportion of seq have to be mismatch after the cuttoffLength
		$minPostCutoffLenASeqPrptn = 0.7; #---proportion of seq have to be A after the cuttoffLength
		$minMisMatchNumPostCutoffLen = 2;#----specify the minimum number of mismatches after cutoffLen
	}

	#---open files for for category of output
	my $NMCutoff = $minMisMatchNumPostCutoffLen;
	my @sRNACtgryAry = ("atLeastNM$NMCutoff.withTail.allSize", "atLeastNM$NMCutoff.noTail.26To28nt", "atLeastNM$NMCutoff.noTail.atLeast29nt","atLeastNM$NMCutoff.noTail.18To25nt", "lesserNM$NMCutoff.26To28nt", "lesserNM$NMCutoff.atLeast29nt","lesserNM$NMCutoff.18To25nt", "lesserNM$NMCutoff.shorter18nt", "atLeastNM$NMCutoff.noTail.shorter18nt");
	foreach my $sRNACtrgy (@sRNACtgryAry) {
		my $outBamPath = "$outDir/$sRNACtrgy.$misMatchMode.$bamPathSplt[-1]";
		$outBamPath =~ s/sam$/bam/;
		$outBamPathHsh{$sRNACtrgy} = $outBamPath;
		$sRNACtgryCountHsh_ref->{$sRNACtrgy} = 0;
		open $samFileHandleHsh{$sRNACtrgy}, "| samtools view -S -b - >$outBamPath 2>>$outDir/log/samtool.log.txt" if ($countMismatchOnly ne 'yes'); 
	}
	my $sRNACtgryAfterSorting = '';
	my %tmpSeqHsh;

	#---print header
	if ($countMismatchOnly eq 'no') {
		open (BAMHEADER, " samtools view -H $bamPath | ");
		while (<BAMHEADER>) {
			chomp; 
			foreach my $sRNACtrgy (@sRNACtgryAry) {
				print {$samFileHandleHsh{$sRNACtrgy}} $_."\n";
			}
		}
		close BAMHEADER;
	}
	
	while (my $theLine = <INBAM>) {

		# HTt83343643	0	EhR1rDNA	11053	255	22M611N16M	*	0	0	ATTCCCACTGTCCCTATCTGCATTTCAAGCAGAATTGA	7:767)8777D=DDDDADCD=@6@?A?AAA@B@@?A>A	XA:i:0	XS:A:+
		
		#---report the progress
		chomp $theLine;
		
		#---check progress
		$lineProc++; $progCount++;
		if ($progCount >= $intervalSize) {
			($progCount, $intervalStart) = reportProgress($progCount, $lineProc, $intervalSize, $fileTotalLineNum, $intervalStart);
		}

		#----ad hoc cmd to stop the 
		#last if $lineProc >= 1000000;
		
		#---split the line if not header
		my (undef, $SAMFlag, $curntCntg, $readStart, undef, $cigarStr, undef, undef, undef, $alignSeq, undef) = split /\t/, $theLine;
		my $readLength = length $alignSeq;

		my $NH = my $NM = my $MD = -1;
		$NM = $1 if ($theLine =~ m/NM\:i\:(\S+)/);
		$NH = $1 if ($theLine =~ m/NH\:i\:(\S+)/);
		$MD = $1 if ($theLine =~ m/MD\:Z\:(\S+)/);
		
		if ($curntCntg eq "*") {
			print "reached unaligned reads. Assuming the SAM file is sorted. Quitting.\n";
			last;
		}#---unaligned
		
		#---check strand
		my $readStrand;
		
		if ($SAMFlag & 16) {
			$readStrand = '-';
		} else {
			$readStrand = '+';
		}

		if ($hardCodedCriteria eq "yes") {
			#next if $curntCntg eq "EhR1rDNA"; #---rRNA
			#next if ($readLength >35);
			#next if (($readLength <28) or ($readLength >35));
		}

		my $rdSeq = $alignSeq;
		if ($readStrand eq "-") {
			$rdSeq = reverse $rdSeq;
			$rdSeq =~ tr/ACGTacgt/TGCAtgca/;
		}

		#---check redundant read
		%tmpSeqHsh = () if (not exists ${$cntgPosHsh{$curntCntg}}{$readStart}); #---clear the tmpSeqHsh if the currnt alignment is at new position

		if (exists $tmpSeqHsh{$rdSeq}) {#---redundant if exists $tmpSeqHsh{$rdSeq}

			#---record the redundant position
			my $read5End = $readStart;
			$read5End = $readStart+$readLength-1 if $readStrand eq '-';
			my $locationTag = $curntCntg.":".$read5End.":".$readStrand;
			$redundantPosHsh{$locationTag}++;

			#--skip the read if removeRedundant
			next if $removeRedundant eq 'yes';
		}

		$tmpSeqHsh{$rdSeq}++;
		${$cntgPosHsh{$curntCntg}}{$readStart}++;
		
		#---set the cutoff length
		my $cutoffLen;
		if ($setCutoffLen eq "freeLen") {
			$cutoffLen = $readLength - 5; #---cutoff Length at 5 nt upstream
		} elsif ($setCutoffLen eq "post26ntATail") {
			$cutoffLen = 26; #---cutoff Length at 26nt 
		}
		
		#---categorize the reads
		my $smallRNAMismatchTail = "no";
		if ($NM >= $NMCutoff) {#----with mismatch
			
			if ($readLength >= $cutoffLen + $minMisMatchNumPostCutoffLen) {
				my $postCutoffLen = $readLength - $cutoffLen;
				my $postCutoffSeq = substr $rdSeq, $cutoffLen, $postCutoffLen;
				my $postCutoffLenASeqNum = $rdSeq =~ tr/A//;
				my $postCutoffLenASeqPrptn = $postCutoffLenASeqNum/$postCutoffLen;
			
				if ($postCutoffLenASeqPrptn >= $minPostCutoffLenASeqPrptn) {
					#actual read 	AGAAATTAGTAGAAAAGAAAAAAAA   MD:Z:0T18G4T0   NM:i:3  NH:i:1
					#genome seq 	TGAAATTAGTAGAAAAGAAGAAAAT	
					my @MDSeqSplt = split /A|T|G|C|N/, $MD; #---(0, 18, 4, 0);
					my @MDNumSplt = split /\d+/, $MD; #---(,T, G, T,);
					my %mismtchPosHsh;
					my $misMatchNumPostCutoffLen = 0;
					my $curntPos = 0;
					my $mismatchNum = @MDNumSplt;
					my @alignSeqSplt = split //, $alignSeq;
					#print $MD."\n";#---debug
					for my $index (0..($mismatchNum-2)) {
						my $numResidue = $MDSeqSplt[$index];
						my $genomeResSeq = $MDNumSplt[$index+1];
						$curntPos += $numResidue;
						my $readResSeq = $alignSeqSplt[$curntPos];
						my $posOnRead = $curntPos+1;
						if ($readStrand eq "-") {
							$readResSeq =~ tr/ACGTacgt/TGCAtgca/;
							$genomeResSeq =~ tr/ACGTacgt/TGCAtgca/;
							$posOnRead = $readLength - $curntPos;
						}
						$misMatchNumPostCutoffLen ++ if $posOnRead > $cutoffLen;
						@{$mismtchPosHsh{$posOnRead}} = ($genomeResSeq, $readResSeq); #--- 0->(T, A), 19->(G, A), 24->(T, A);
						#print $readStrand."\t".$readLength."\t".$posOnRead."\t".$genomeResSeq."\t".$readResSeq."\n";#---debug
						$curntPos++;
					}
				
					if ($misMatchNumPostCutoffLen >= $minMisMatchNumPostCutoffLen) {
						my $misMatchAsPostCutoffLen = 0;
						my @misMatchAPosAry;
						foreach my $posOnRead (sort {$b <=> $a} keys %mismtchPosHsh) {
							my ($genomeResSeq, $readResSeq) = @{$mismtchPosHsh{$posOnRead}};
							if (($readResSeq eq "A") or ($readResSeq eq "a")) {
								push @misMatchAPosAry, $posOnRead;
								$misMatchAsPostCutoffLen++;
							}
						}
					
						my $misMatchAsPostCutoffLenPrptn = $misMatchAsPostCutoffLen/$misMatchNumPostCutoffLen;
						if ($misMatchAsPostCutoffLenPrptn >= $minPostCutoffLenAMismatchPrptn) {
							#print $readStrand."\n";
							$smallRNAMismatchTail = 'yes';
							foreach my $posOnRead (sort {$b <=> $a} keys %mismtchPosHsh) {
								my ($genomeResSeq, $readResSeq) = @{$mismtchPosHsh{$posOnRead}};
								if ($posOnRead > $cutoffLen) {
									${${$misMatchCountHsh{"postCutoffLen"}}{$readResSeq}}{$genomeResSeq}++;
								} else {
									${${$misMatchCountHsh{"preCutoffLen"}}{$readResSeq}}{$genomeResSeq}++;
								}
							}
						}
					}
				}
			}
			
			if ($smallRNAMismatchTail eq 'yes') {
				$sRNACtgryAfterSorting = "atLeastNM$NMCutoff.withTail.allSize";
			} else {
				if ($readLength >= 26 and $readLength <= 28) {
					$sRNACtgryAfterSorting = "atLeastNM$NMCutoff.noTail.26To28nt";
				} elsif ($readLength >= 18 and $readLength <= 25) {
					$sRNACtgryAfterSorting = "atLeastNM$NMCutoff.noTail.18To25nt";
				} elsif ($readLength >= 29) {
					$sRNACtgryAfterSorting = "atLeastNM$NMCutoff.noTail.atLeast29nt";
				} else {
					$sRNACtgryAfterSorting = "atLeastNM$NMCutoff.noTail.shorter18nt";
				}
			}

		} else {#---end of if (($NM >= $NMCutoff) and ($readLength > $cutoffLen))
			
			if ($readLength >= 26 and $readLength <= 28) {
				$sRNACtgryAfterSorting = "lesserNM$NMCutoff.26To28nt";
			} elsif ($readLength >= 18 and $readLength <= 25) {
				$sRNACtgryAfterSorting = "lesserNM$NMCutoff.18To25nt";
			} elsif ($readLength >= 29){
				$sRNACtgryAfterSorting = "lesserNM$NMCutoff.atLeast29nt";
			} else {
				$sRNACtgryAfterSorting = "lesserNM$NMCutoff.shorter18nt";
			}
		}
		print {$samFileHandleHsh{$sRNACtgryAfterSorting}} $theLine."\n" if ($countMismatchOnly eq 'no');
		$sRNACtgryCountHsh_ref->{$sRNACtgryAfterSorting}++;
		$sRNACtgryCountHsh_ref->{'total'}++;
	}

	print  "\n";
	
	close INBAM;
	
	if ($countMismatchOnly eq 'no') {
		foreach my $sRNACtrgy (@sRNACtgryAry) {
			close $samFileHandleHsh{$sRNACtrgy}; 
		}
	}
	
	#----print the sorted counts on screen
	open (CATCOUNT, ">", "$outDir/log/categoryCount.log.txt");
	foreach my $sRNACtrgy (sort keys %{$sRNACtgryCountHsh_ref}) {
		my $pct = sprintf "%.2f", 100*$sRNACtgryCountHsh_ref->{$sRNACtrgy}/$sRNACtgryCountHsh_ref->{'total'};
		my $outline = sprintf ("%30s : %-12s : %-10s", $sRNACtrgy, $sRNACtgryCountHsh_ref->{$sRNACtrgy}, $pct."%");
		print $outline."\n";
		print CATCOUNT $outline."\n";
		delete $outBamPathHsh{$sRNACtrgy} if $sRNACtgryCountHsh_ref->{$sRNACtrgy} == 0;
	}
	close CATCOUNT;
	
	return \%outBamPathHsh, \%misMatchCountHsh, \%redundantPosHsh;
}
########################################################################## checkFileSizeAndDefineIntervalSize
sub checkFileSizeAndDefineIntervalSize {
    
    my $fileToCheckPath = $_[0];
    my $linesToSample = $_[1];
	
	#---make sure $linesToSample is a non-zero number, if not set to 10000
	my $linesToSampleInt = int $linesToSample;
	$linesToSampleInt = 100000 if (($linesToSampleInt != $linesToSample) or ($linesToSampleInt == 0));

    #---get the filename from the path
    my @fileToCheckPathSplt = split /\//, $fileToCheckPath;
    my $fileToCheckName = $fileToCheckPathSplt[-1];
    
    print "Estimating the number of lines in $fileToCheckName.\n";
	
	#---estimate the number of lines in the file
	open (INFILE, $fileToCheckPath) || die "Can't open $fileToCheckPath.\n";close INFILE;
	my $tmpFilePath = $fileToCheckPath."_tmp.txt";
	system "tail -$linesToSampleInt $fileToCheckPath >$tmpFilePath";
	my $fileToCheckSize = -s "$fileToCheckPath";
   	my $tmpFileSize = -s "$tmpFilePath";
	system "rm $tmpFilePath";
   	my $fileToCheckSizeTotalLineNum = int (($fileToCheckSize/$tmpFileSize)*$linesToSampleInt);
   	print "Estimated to have ".$fileToCheckSizeTotalLineNum." lines in $fileToCheckName.\n";

	my $intervalSize = int ($fileToCheckSizeTotalLineNum/100); #---define as 
	$intervalSize = 1000000 if  ($intervalSize > 1000000);
	
	return ($fileToCheckSizeTotalLineNum, $intervalSize);

}
########################################################################## reportProgress
sub reportProgress {

	my $progCount = $_[0];
	my $lineProc = $_[1];
	my $intervalSize = $_[2];
	my $fileTotalLineNum = $_[3];
	my $intervalStart = $_[4];

	$progCount=0;
	my $intervalEnd = time();
	my $timeElapsed = $intervalEnd - $intervalStart;
	$timeElapsed = sprintf ("%.2f", $timeElapsed);
	my $estimatedEnd = (($fileTotalLineNum - $lineProc)*$timeElapsed)/$intervalSize;
	$estimatedEnd = sprintf ("%.2f", $estimatedEnd/60);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	updateProgressBar("[$runTime]; Interval=$timeElapsed sec; EstimatedEnd=$estimatedEnd min", $lineProc, $fileTotalLineNum, 50, 5);
	#print "$lineProc lines processed. Last $intervalSize lines:".$timeElapsed." sec. Estimated end: ".$estimatedEnd." mins.\r";
	$intervalStart = time();
	
	return ($progCount, $intervalStart);
		
}
########################################################################## updateProgressBar
sub updateProgressBar {
	
	my $strToPrint = $_[0];
	my $progressCount = $_[1];
	my $totalCount = $_[2];
	my $scaleMax = $_[3];
	my $extraWhiteSpaceNum = $_[4]; #---erase the longer infos during the progress
	
	my $progressPct = int (($progressCount/$totalCount)*$scaleMax);

	my $progressBar = "|";
	for my $i (1..$progressPct) {$progressBar .= ">";}
	for my $i (($progressPct+1)..$scaleMax) {$progressBar .= " ";}
	$progressBar .= "|";

	my $extraWhiteSpaceStr = "";
	for my $i (1..$extraWhiteSpaceNum) {$extraWhiteSpaceStr .= " ";}
	
	print $progressBar.$strToPrint.$extraWhiteSpaceStr."\r";

}
########################################################################## printProgressScale
sub printProgressScale {

	my $strToPrint = $_[0];
	my $scaleMax = $_[1];

	my $scaleSpace = "|";
	for my $i (1..$scaleMax) {$scaleSpace .= "-";}
	$scaleSpace .= "|100%";
	
	print $strToPrint."\n";
	print $scaleSpace."\n";
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## runSamtoolIndexBAM
sub runSamtoolIndexBAMAndPileupCounter {
	
	my ($outBamPathHsh_ref, $pileupCounterPath, $GFFPath, $fastaPath) = @_;
	my %outBamPathHsh = %{$outBamPathHsh_ref};
	
	#---run samtools
	print "Indexing the bam files\n";
	foreach my $sRNACtrgy (keys %outBamPathHsh) {
		my $cmd = "samtools index $outBamPathHsh{$sRNACtrgy}";
		system (qq|$cmd|);
	}
	
	#---run samtools
	my %cmdHsh;
	print "Indexing the bam files\n";
	foreach my $sRNACtrgy (keys %outBamPathHsh) {
		my $cmd = "perl $pileupCounterPath --pileupPath=no --gffPath=$GFFPath --printReadSeq=yes --refFastaPath=$fastaPath --covCutoff=1 --breakdownSam=no --cntgFiltr=all --bamPath=$outBamPathHsh{$sRNACtrgy} --outDir=$outDir/pileupCounter/$sRNACtrgy/";
		$cmdHsh{$sRNACtrgy} = $cmd;
		print "Issued pileupCounter on $sRNACtrgy\n";
		system (qq|$cmd 1>$outDir/log/pileupCounterRunlog/$sRNACtrgy.log.txt &|);
	}

	print "\n";
	
	my $pileupCounterFinish = 'no';
	while ($pileupCounterFinish eq 'no') {
		my @sRNACtrgyRunningAry = ();
		foreach my $sRNACtrgy (keys %cmdHsh) {
			my $cmd = $cmdHsh{$sRNACtrgy};
			my $stdout = `ps -ef | grep -e '$cmd' | grep -v grep`;
			push @sRNACtrgyRunningAry, $sRNACtrgy if $stdout =~ m/$cmd/;
		}
		sleep 10;
		my $sRNACtrgyStr = join ", " , @sRNACtrgyRunningAry;
		my $outLine = sprintf "%-200s", "pileupCounter for $sRNACtrgyStr are still running";
		print $outLine."\r";
		$pileupCounterFinish = 'yes' if @sRNACtrgyRunningAry == 0;
	}

	print "\n";
	print "pileupCounter finished\n";

}
########################################################################## printMisMatchCount
sub printMisMatchCount {
	
	#---printMisMatchCount($misMatchCountHsh_ref, $outDir);
	my %misMatchCountHsh = %{$_[0]};
	my $misMatchMode = $_[1];
	my $outDir = $_[2];
	
	open (MISMATCHLOG, ">$outDir/log/$misMatchMode.tail.mismtach.log.txt");
	print MISMATCHLOG join "", ((join "\t", ('preOrPostCutoffLen', 'readResSeq', 'genomeResSeq', 'count')), "\n");
	foreach my $preOrPostCutoffLen (keys %misMatchCountHsh ) {
		foreach my $readResSeq (sort {$a cmp $b} keys %{$misMatchCountHsh{$preOrPostCutoffLen}}) {
			foreach my $genomeResSeq (sort {$a cmp $b} keys %{${$misMatchCountHsh{$preOrPostCutoffLen}}{$readResSeq}}) {
				my $count = ${${$misMatchCountHsh{$preOrPostCutoffLen}}{$readResSeq}}{$genomeResSeq};
				print MISMATCHLOG join "", ((join "\t", ($preOrPostCutoffLen, $readResSeq, $genomeResSeq, $count)), "\n");
			}
		}
	}
	close MISMATCHLOG;
}
########################################################################## printMisMatchCount
sub printRedundantPos {
	
	#---printMisMatchCount($misMatchCountHsh_ref, $outDir);
	my %redundantPosHsh = %{$_[0]};
	my $outDir = $_[1];
	
	open (RDNDNTLOG, ">$outDir/log/redundant.readEnd5Pos.log.txt");
	print RDNDNTLOG join "", ((join "\t", ('locationTag', 'count')), "\n");
	foreach my $locationTag (sort {$redundantPosHsh{$b} <=> $redundantPosHsh{$a}} keys %redundantPosHsh ) {
		my $count = $redundantPosHsh{$locationTag};
		print RDNDNTLOG join "", ((join "\t", ($locationTag, $count)), "\n");
	}
	close RDNDNTLOG;
}
########################################################################## checkBAMAlignNum
sub checkBAMAlignNum {
	
	my $bamPathToCheck = $_[0];
	
	print "Checking the number of aligments in the bam file, please be patient, it takes a while.\n";
	
	my $flagStatOut = `samtools flagstat $bamPathToCheck`;
	
	my @firstLineSplt = split / /, $flagStatOut;
	my $fileToCheckSizeTotalLineNum = $firstLineSplt[0];
	
	die "Problem in check bam file size\n" if ($fileToCheckSizeTotalLineNum < 1);
	
	my $intervalSize = int ($fileToCheckSizeTotalLineNum/100); #---define as 
	$intervalSize = 1000000 if  ($intervalSize > 1000000);

	return ($fileToCheckSizeTotalLineNum, $intervalSize);
}
