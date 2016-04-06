#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;

use MyConfig;
use Bed;
use ScanStatistics;
use Data::Dumper;

my $cache = 'cache/';
my $exonTagCountBedFile = "$cache/tag.count.exon.bed";
my $fin;
my $fout;

my $verbose = 1;

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

print Dumper(%geneTagCountHash);
my $geneTagCountTotal = 0;
map {$geneTagCountTotal += $geneTagCountHash{$_}->{'count'}} keys %geneTagCountHash;

Carp::croak "no tags??" if $geneTagCountTotal <= 0;

print "total number of tags overlapping with specified reions: $geneTagCountTotal\n" if $verbose;
