#!/usr/bin/env perl

#date:   2025-03-11
#writer: Xiangjian Gou (QQ = 862137261)

use strict;
use warnings;

my ($samtoolsDir, $gatkDir, $outputDir) = @ARGV;

open my $samtools_snp_fh, '<', "$samtoolsDir/samtools.SNP.vcf" or die "Error: cannot open samtools.SNP.vcf: $!";
open my $gatk_snp_fh, '<', "$gatkDir/gatk.SNP.vcf" or die "Error: cannot open gatk.SNP.vcf: $!";
open my $samtools_indel_fh, '<', "$samtoolsDir/samtools.INDEL.vcf" or die "Error: cannot open samtools.INDEL.vcf: $!";
open my $gatk_indel_fh, '<', "$gatkDir/gatk.INDEL.vcf" or die "Error: cannot open gatk.INDEL.vcf: $!";

foreach my $variant ( [$samtools_indel_fh, $gatk_indel_fh], [$samtools_snp_fh, $gatk_snp_fh] ) {
    my %name_info;
    foreach my $i (0, 1) {
        my $tmp_fh = $variant->[$i]; #cannot use $variant->[$i] in <>
        while (<$tmp_fh>) {
            if (/\A#/) {
                $name_info{gatk_note} .= $_ if $i == 1;
                next;
            }
            my $name = join ':', (split)[0, 1, 3, 4];
            $name_info{$name}[0]++;
            $name_info{$name}[1] = $_; #GATK will cover samtools
        }
    }

    #output the result
    open my $out_fh, '>', $variant->[0] eq $samtools_indel_fh ? "$outputDir/INDEL.vcf" : "$outputDir/SNP.vcf";
    print $out_fh $name_info{gatk_note};
    delete $name_info{gatk_note};
    my %tmp;
    foreach my $name (keys %name_info) {
        next if $name_info{$name}[0] != 2;
        my ($chr, $loci) = (split /:/, $name)[0, 1];
        $tmp{$chr}{$loci} = $name_info{$name}[1];
    }
    undef %name_info;
    foreach my $chr (sort keys %tmp) {
        foreach my $loci (sort {$a <=> $b} keys %{$tmp{$chr}}) {
            print $out_fh $tmp{$chr}{$loci};
        }
    }
    close $out_fh;
}

close $samtools_snp_fh;
close $gatk_snp_fh;
close $samtools_indel_fh;
close $gatk_indel_fh;
