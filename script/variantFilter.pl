#!/usr/bin/env perl

#date:   2025-03-11
#writer: Xiangjian Gou (QQ = 862137261)

use strict;
use warnings;

$| = 1;

my ($inputFile, $mafThreshold, $mrThreshold, $outputFilePrefix) = @ARGV;

open my $iHMP, '<', $inputFile or die "Error: cannot open file '$inputFile': $!";
open my $oHMP, '>', "$outputFilePrefix-MAF$mafThreshold-MR$mrThreshold.hmp";
print $oHMP scalar <$iHMP>;

while (<$iHMP>) {
    my @row = split;
    my ($snpID, $allele) = @row[0, 1];
    next if length($allele) != 3; #only keep bi-allelic SNPs

    #count the number of each genotype
    my %baseCount;
    my $hetCount = 0;
    my $nCount = 0;
    my @snp = @row[11 .. $#row];
    foreach my $base (@snp) {
        if ($base eq 'A' or $base eq 'T' or $base eq 'G' or $base eq 'C' or $base eq '-') {
            $baseCount{$base}++;
        }
        elsif ($base eq 'N') {
            $nCount++;
        }
        else {
            $hetCount++;
        }
    }

    #judge and output
    my ($maxCount, $minCount);
    if (keys %baseCount == 0) {
        $maxCount = 0;
        $minCount = 0;
    }
    elsif (keys %baseCount == 1) {
        my $maxBase = (sort keys %baseCount)[0];
        $maxCount = $baseCount{$maxBase};
        $minCount = 0;
    }
    elsif (keys %baseCount == 2) {
        my ($maxBase, $minBase) = sort { $baseCount{$b} <=> $baseCount{$a} } keys %baseCount;
        $maxCount = $baseCount{$maxBase};
        $minCount = $baseCount{$minBase};
    }
    else {
        die "Error: bug1 in hapmap file ($snpID) !\n";
    }
    next if $maxCount+$minCount+$hetCount == 0;

    #filter
    my $MAF = sprintf "%.3f", ($minCount+0.5*$hetCount)/($maxCount+$minCount+$hetCount); #MAF: exclude N
    my $MR  = sprintf "%.3f", $nCount/@snp; #MR: only consider N
    print $oHMP $_ if $MAF >= $mafThreshold and $MR <= $mrThreshold;
}

close $iHMP;
close $oHMP;
