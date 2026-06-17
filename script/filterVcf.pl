#!/usr/bin/env perl

#date:   2025-03-11
#writer: Xiangjian Gou (QQ = 862137261)

use strict;
use warnings;

$| = 1;

my ($vcfFile, $outputFile1, $outputFile2) = @ARGV;

open my $iVCF, '<', $vcfFile or die "Error: cannot open file '$vcfFile': $!";
open my $oVCF1, '>', $outputFile1; #ref base is not het
open my $oVCF2, '>', $outputFile2; #ref base is het
while (<$iVCF>) {
    if (/\A#/) {
        print $oVCF1 $_;
    }
    else {
        my $ref = (split)[3];
        $ref =~ s/[AGTCN]//g;
        if (! $ref) {
            print $oVCF1 $_;
        }
        else {
            print $oVCF2 $_;
        }
    }
}
close $iVCF;
close $oVCF1;
close $oVCF2;
