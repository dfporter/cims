#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;

use MyConfig;
use Bed;
use ScanStatistics;


my $prog = basename ($0);
my $pvalueThreshold = 0.01;
my $maxGap = 10;
my $cache = getDefaultCache ($prog);
my $keepCache = 0;
my $verbose = 0;
my $big = 0;
my $separateStrand = 0;
my $useExpr = 0;
my $multiTestCorrection = 0;

my @ARGV0 = @ARGV;

GetOptions (
	'big'=>\$big,
	'p:f'=>\$pvalueThreshold,
	'ss'=>\$separateStrand,
	'gap:i'=>\$maxGap,
	'use-expr'=>\$useExpr,
	'multi-test'=>\$multiTestCorrection,
	'c|cache:s'=>\$cache,
	'keep-cache'=>\$keepCache,
	'v'=>\$verbose);


if (@ARGV != 3)
{
    print "detecting significant peaks from CLIP data\n";
    print "Usage: $prog [options] <gene.bed> <tag.bed> <out.bed>\n";
    #print " -e          : expression values indicated in the gene bed file\n";
    print " -big        : big file\n";
	print " -ss         : separate the two strands\n";
    print " --use-expr  : use expression levels given in the score column in the gene bed file for normalization\n";
	print " -p   [float]: threshold of pvalue to call peak ($pvalueThreshold)\n";    
	print " --multi-test: do bonforroni multiple test correction\n";
    print " -gap   [int]: cluster peaks closer than the gap ($maxGap)\n";
    print " -c     [dir]: cache dir\n";
	print " --keep-cache: keep cache when the job done\n";
    print " -v          : verbose\n";
    exit (0);
}

print "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;

my ($geneBedFile, $tagBedFile, $outBedFile) = @ARGV;

my $cmdDir = dirname ($0);

my $bigFlag = $big ? "-big" : "";
my $verboseFlag = $verbose ? "-v" : "";

system ("mkdir $cache") unless -d $cache;
print "estimating average tag size ...\n" if $verbose;

my $tagSize =  `awk 'BEGIN{s=0;n=0;} {s=s+\$3-\$2; n=n+1} END {print s/n}' $tagBedFile`;
chomp $tagSize;
print "average tag size = $tagSize\n" if $verbose;

print "get exons in genes ... \n" if $verbose;
my $ts2geneFile = "$cache/ts2gene.txt";
my $cleanGeneBedFile = "$cache/gene.clean.bed";

my $cmd = "awk '{print \$4\"\\t\"\$4}' $geneBedFile | sort | uniq > $ts2geneFile";
system ($cmd);

$cmd = "perl $cmdDir/combineTranscripts.pl $verboseFlag  $geneBedFile $ts2geneFile $cleanGeneBedFile";
system ($cmd);

my $exonBedFile = "$cache/exon.bed";
$cmd = "perl $cmdDir/gene2ExonIntron.pl $verboseFlag -nid -oe $exonBedFile $cleanGeneBedFile";
system ($cmd);


print "count tag number for each exon/gene ...\n" if $verbose;
my $exonTagCountBedFile = "$cache/tag.count.exon.bed";

my $ssFlag = $separateStrand ? "-ss" : "";
$cmd = "perl $cmdDir/tag2profile.pl $bigFlag $verboseFlag -region $exonBedFile $ssFlag -of bed $tagBedFile $exonTagCountBedFile";
print $cmd, "\n";
system ($cmd);

my $fin;
my $fout;

my %geneTagCountHash;

open ($fin, "<$exonTagCountBedFile") || Carp::croak "cannot open file $exonTagCountBedFile\n";
while (my $line = <$fin>)
{
    chomp $line;
    my $r = lineToBed ($line);

    my $geneId = $r->{"name"};
    $geneTagCountHash{$geneId}->{'size'} += ($r->{'chromEnd'} - $r->{'chromStart'} + 1);
    $geneTagCountHash{$geneId}->{'count'} += $r->{'score'};
}
close ($fin);

my $geneTagCountTotal = 0;
map {$geneTagCountTotal += $geneTagCountHash{$_}->{'count'}} keys %geneTagCountHash;

Carp::croak "no tags??" if $geneTagCountTotal <= 0;

print "total number of tags overlapping with specified reions: $geneTagCountTotal\n" if $verbose;

my $geneExpressionTotal = 0;
if ($useExpr)
{
	print "reading gene expression levels from $geneBedFile ...\n" if $verbose;
	open ($fin, "<$geneBedFile") || Carp::croak "cannot open file $geneBedFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^track/;
		next if $line =~/^#/;

		my @cols = split (/\s+/, $line);
		Carp::croak "no score column in $geneBedFile ...\n" if @cols < 5;
		
		my $geneId = $cols[3];
		my $expr = $cols[4];
		
		Carp::croak "expression value for gene $geneId is negative\n" if $expr < 0;
		$geneTagCountHash{$geneId}->{'expr'} = $expr if exists $geneTagCountHash{$geneId};
	}
	close ($fin);

	map {$geneExpressionTotal += $geneTagCountHash{$_}->{'expr'}} keys %geneTagCountHash;
	Carp::croak "all genes have no expression??\n" if $geneExpressionTotal <= 0;

	print "total gene expression level: $geneExpressionTotal\n" if $verbose;
}

print "calculating expected peak height for each gene ...\n" if $verbose;


#my $geneTagCountFile = "$cache/tag.count.gene.txt";
#open ($fout, ">$geneTagCountFile") || Carp::croak "cannot open file $geneTagCountFile to write\n";

#if ($useExpr)
#{
#	print $fout "#", join ("\t", "gene id", "size", "tag number", "tag density", "expression", "expexpted PH"), "\n";
#}
#else
#{
#	print $fout "#", join ("\t", "gene id", "size", "tag number", "tag density", "expexpted PH"), "\n";
#}

foreach my $geneId (sort keys %geneTagCountHash)
{
	my $g = $geneTagCountHash{$geneId};
    $g->{'density'} = $g->{'count'} / $g->{'size'};
	$g->{'PH0'} = $useExpr ? ($geneTagCountTotal * $g->{'expr'} / $geneExpressionTotal / $g->{'size'} * $tagSize) : ($g->{'density'} * $tagSize);

#	if ($useExpr)
#	{
#		print $fout join ("\t", $geneId, $g->{'size'}, $g->{'count'}, $g->{'density'}, $g->{'expr'}, $g->{'PH0'}), "\n";
#	}
#	else
#	{
#		print $fout join ("\t", $geneId, $g->{'size'}, $g->{'count'}, $g->{'density'}, $g->{'PH0'}), "\n";
#	}
}
#close ($fout);

#my $effectiveGeneNum = `awk '{if(\$3>0) {print \$0}}' $geneTagCountFile | wc -l`;
#$effectiveGeneNum =~/\s*(\d+)\s*/;
#$effectiveGeneNum = $1;

my $effectiveGeneNum = 0;

if ($useExpr)
{
	map {$effectiveGeneNum++ if $geneTagCountHash{$_}->{'expr'} > 0} keys %geneTagCountHash;
}
else
{
	map {$effectiveGeneNum++ if $geneTagCountHash{$_}->{'count'} > 0} keys %geneTagCountHash;
}

print "effective gene number: $effectiveGeneNum\n" if $verbose;

print "generate exact tag profile ...\n" if $verbose;

my $tagExactCountBedFile = "$cache/tag.count.exact.bed";

$cmd = "perl $cmdDir/tag2profile.pl $verboseFlag $bigFlag -exact $ssFlag -of bed $tagBedFile $cache/tmp.bed";
system ($cmd);

$cmd = "awk '{if (\$5>1) {print \$0}}' $cache/tmp.bed > $tagExactCountBedFile";
system ($cmd);

$cmd = "wc -l $tagExactCountBedFile";
my $n = `$cmd`;
$n=~/^(\d+)/;
$n = $1;

Carp::croak "no peaks found\n" if $n == 0;

print "$n positions with more than one tag detected\n" if $verbose;


print "match peaks with genes ...\n" if $verbose;
my $peak2geneBedFile = "$cache/peak2gene.bed";
$cmd = "perl $cmdDir/tagoverlap.pl $bigFlag $verboseFlag -d \"@@\" --keep-score -region $cleanGeneBedFile $ssFlag $tagExactCountBedFile $peak2geneBedFile";
system ($cmd);


print "calculate pvalues ...\n" if $verbose;

open ($fin, "<$peak2geneBedFile") || Carp::croak "cannot open file $peak2geneBedFile to read\n";
open ($fout, ">$cache/tmp.bed") || Carp::croak "cannot open file $cache/tmp.bed to write\n";

my %resultHash;
my $i = 0;

while (my $line = <$fin>)
{
    chomp $line;

    print "$i ...\n" if $verbose && $i % 5000 == 0;
    $i++;

    my $peak = lineToBed ($line);
    my $name = $peak->{'name'};
    my $peakHeight = $peak->{'score'};

    my ($peakId, $geneId) = split ("@@", $name);
    next unless exists $geneTagCountHash{$geneId};

    my $expectedPeakHeight = $geneTagCountHash{$geneId}->{'PH0'}; #expected PH
    next unless $expectedPeakHeight > 0; #otherwise no tag on the exonic region of the gene
    next unless $peakHeight > $expectedPeakHeight;#otherwise it cannot be significant

    my $geneSize = $geneTagCountHash{$geneId}->{'size'};

    my $pvalue = exists $resultHash{$geneId}{$peakHeight} ? 
    	$resultHash{$geneId}{$peakHeight} : calcScanStatistic ($peakHeight, $expectedPeakHeight, $geneSize / $tagSize);

    $resultHash{$geneId}{$peakHeight} = $pvalue unless exists $resultHash{$geneId}{$peakHeight};

    $pvalue *= $effectiveGeneNum if $multiTestCorrection;
    $pvalue = 1e-100 if $pvalue == 0;

    $peak->{'name'} = $peakId . "[gene=$geneId][PH=$peakHeight][PH0=" . 
    	sprintf ("%.2f", $expectedPeakHeight) . "][P=" .sprintf ("%.2e", $pvalue) . "]";
    $peak->{'score'} = $peakHeight; #-log ($pvalue) / log(10);
    
    print $fout bedToLine ($peak), "\n" if $pvalue <= $pvalueThreshold;
    
}

close ($fin);
close ($fout);

print "clustering peaks ...\n" if $verbose;
$cmd = "perl $cmdDir/bedUniq.pl $ssFlag -c maxScore -maxgap $maxGap $verboseFlag $cache/tmp.bed $outBedFile";
system ($cmd);

system ("rm -rf $cache") unless $keepCache;
