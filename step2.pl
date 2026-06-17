#!/usr/bin/env perl

#date    : 2026-06-17
#writer1 : Xiangjian Gou (xjgou@mail.hzau.edu.cn)
#writer2 : Haoran Shi (sicau3339@outlook.com)

#load modules
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/getcwd abs_path/;

#record version information
my $VERSION = 'VCP, step2 v1.3 (2026-06-17)';

#get the home directory of VCP
my $vcp_path = abs_path $0;
$vcp_path =~ s/\\/\//g;
$vcp_path =~ s/\/[^\/]+\z//;
my $softwarePath = "$vcp_path/required_software.txt";
my $intersectionScript = "$vcp_path/script/getIntersectionVariants.pl";
my $filterScript = "$vcp_path/script/filterVcf.pl";

#get current working directory
my $cwd = getcwd;

#set default options
my $germplasmDir = "/public/home/xjgou/raw_reads";
my $refGenomeFa = "/public/home/xjgou/genome/genome.fa";
my $outputDir = "output";
my $jumpSamtools = '';
my $jumpGatk = '';
my $jumpIntersection = '';
my $queue = "smp";
my $thread = 10;
my $memory = 100;
my $version;
my $help;

#get options from command line
GetOptions(
    'germplasmDir=s'        => \$germplasmDir,
    'refGenomeFa=s'         => \$refGenomeFa,
    'outputDir=s'           => \$outputDir,
    'jumpSamtools|js+'      => \$jumpSamtools,
    'jumpGatk|jg+'          => \$jumpGatk,
    'jumpIntersection|ji+'  => \$jumpIntersection,
    'queue=s'               => \$queue,
    'thread=i'              => \$thread,
    'memory=i'              => \$memory,
    'version+'              => \$version,
    'help+'                 => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Name: VCP - Variant Calling Pipeline (step2)

Function: known sites construction for subsequent Base Quality Score Recalibrator

Usage: perl step2.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for path:
  -g  | -germplasmDir <STRING> : the directory of storing all germplasms. (default: /public/home/xjgou/raw_reads)
  -r  | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/genome/genome.fa)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)

  #Options for step:
  -js | -jumpSamtools          : no execute samtools
  -jg | -jumpGatk              : no execute GATK
  -ji | -jumpIntersection      : no execute intersection

  #Options for resources:
  -q  | -queue  <STRING>       : set the queue to use. (default: smp)
  -t  | -thread    <INT>       : set the number of threads to use. (default: 10)
  -m  | -memory    <INT>       : set the size of memory to use. (default: 100 [100GB])

  #Options for other:
  -v  | -version               : show the version information.
  -h  | -help                  : show the help information.
####################################################################################################

__GUIDE__

#output version and help information
die "$VERSION\n" if $version;
die $usage if $help;

#set whether to execute each step
$jumpSamtools = '#' if $jumpSamtools;
$jumpGatk = '#' if $jumpGatk;
$jumpIntersection = '#' if $jumpIntersection;

#get the all sequence names for subsequent GenomicsDBImport
open my $iFASTA, '<', $refGenomeFa or die "Error: cannot open file '$refGenomeFa': $!";
my @names;
while (<$iFASTA>) {
    if (/\A>/) {
        my ($name) = /\A>(\S+)/;
        push @names, $name;
    }
}
close $iFASTA;
my $chrSets = join " ", map {"-L $_"} @names;

#get the all sorted bam files and gvcf files for subsequent bcftools and GATK, respectively
my @sortedBamFiles;
my @gvcfFiles;
foreach my $dir (glob "$germplasmDir/*") {
    my ($germplasm) = $dir =~ /([^\/]+)\z/;
    push @sortedBamFiles, "$outputDir/SortSamOut/$germplasm.sort.bam";
    push @gvcfFiles, "$outputDir/HaplotypeCallerOut_nonBQSR/$germplasm.g.vcf";
}
my $sortedBamSet = join " ", @sortedBamFiles;
my $gvcfSets = join " ", map {"-V $_"} @gvcfFiles;

#generate a scheduling script for each germplasm, meanwhile, generate a comprehensive scheduling script
my $lsfDir = "step2.lsf";
system "mkdir -p $lsfDir";
open my $oTOTAL, '>', 'step2.sh';
print $oTOTAL "bsub < $lsfDir/1.lsf\n";
close $oTOTAL;
open my $oEACH, '>', "$lsfDir/1.lsf";
lsfInfo($oEACH, $queue, $thread, $memory, $outputDir, $softwarePath, $refGenomeFa, $sortedBamSet, $gvcfSets, $chrSets, $intersectionScript, $filterScript, $jumpSamtools, $jumpGatk, $jumpIntersection, $cwd);
close $oEACH;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $queue, $thread, $memory, $outputDir, $softwarePath, $refGenomeFa, $sortedBamSet, $gvcfSets, $chrSets, $intersectionScript, $filterScript, $jumpSamtools, $jumpGatk, $jumpIntersection, $cwd) = @_;
    my $software = getSoftwarePath($softwarePath);
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -J step2
#BSUB -o step2.out
#BSUB -e step2.err

#create output directory
mkdir -p $outputDir
mkdir -p $outputDir/knownsites

#samtools call
$jumpSamtools mkdir -p $outputDir/knownsites/samtools
$jumpSamtools $software->{BCFTOOLS} mpileup -d 100000 -Ou -f $refGenomeFa $sortedBamSet --threads $thread | $software->{BCFTOOLS} call -m -v -Ob -o $outputDir/knownsites/samtools/samtools.bcf
$jumpSamtools $software->{BCFTOOLS} view $outputDir/knownsites/samtools/samtools.bcf > $outputDir/knownsites/samtools/samtools.vcf
$jumpSamtools perl $filterScript $outputDir/knownsites/samtools/samtools.vcf $outputDir/knownsites/samtools/samtools.refisnothet.vcf $outputDir/knownsites/samtools/samtools.refishet.vcf
$jumpSamtools $software->{GATK} IndexFeatureFile -F $outputDir/knownsites/samtools/samtools.refisnothet.vcf
$jumpSamtools $software->{GATK} VariantFiltration -V $outputDir/knownsites/samtools/samtools.refisnothet.vcf -O $outputDir/knownsites/samtools/samtools.filter.vcf --cluster 4 --window 10 --mask-extension 3 --filter-name "lowMQ" --filter "MQ < 40.0" --filter-name "lowDP" --filter "DP < 8.0" --filter-name "LowQual" --filter "QUAL < 20"
$jumpSamtools $software->{BCFTOOLS} view -f PASS $outputDir/knownsites/samtools/samtools.filter.vcf > $outputDir/knownsites/samtools/tmp
$jumpSamtools mv $outputDir/knownsites/samtools/tmp $outputDir/knownsites/samtools/samtools.filter.vcf
$jumpSamtools rm $outputDir/knownsites/samtools/samtools.filter.vcf.idx
$jumpSamtools $software->{GATK} SelectVariants -O $outputDir/knownsites/samtools/samtools.SNP.vcf --variant $outputDir/knownsites/samtools/samtools.filter.vcf -select-type SNP
$jumpSamtools $software->{GATK} SelectVariants -O $outputDir/knownsites/samtools/samtools.INDEL.vcf --variant $outputDir/knownsites/samtools/samtools.filter.vcf -select-type INDEL

#GATK call
$jumpGatk mkdir -p $outputDir/knownsites/gatk
$jumpGatk $software->{GATK} --java-options "-Xmx${memory}g -Xms${memory}g" GenomicsDBImport --genomicsdb-workspace-path $outputDir/knownsites/gatk/db --batch-size 50 -R $refGenomeFa $gvcfSets $chrSets
$jumpGatk cd $outputDir/knownsites/gatk
$jumpGatk $software->{GATK} --java-options "-Xmx${memory}g -Xms${memory}g" GenotypeGVCFs -R $refGenomeFa -V gendb://db -O gatk.vcf -new-qual -G StandardAnnotation --use-new-qual-calculator
$jumpGatk cd $cwd
$jumpGatk $software->{GATK} VariantFiltration -V $outputDir/knownsites/gatk/gatk.vcf -O $outputDir/knownsites/gatk/gatk.HC.vcf --cluster 4 --window 10 --mask-extension 3 --filter-name "lowMQ" --filter "MQ < 40.0" --filter-name "lowDP" --filter "DP < 8.0" --filter-name "LowQual" --filter "QUAL < 20" --filter-name "lowQD" --filter "QD < 2.0" --filter-name "lowReadPosRankSum" --filter "ReadPosRankSum < -8.0" --filter-name "highFS" --filter "FS > 60.0" --filter-name "lowMQRankSum" --filter "MQRankSum < -12.5"
$jumpGatk $software->{BCFTOOLS} view -f PASS $outputDir/knownsites/gatk/gatk.HC.vcf > $outputDir/knownsites/gatk/tmp
$jumpGatk mv $outputDir/knownsites/gatk/tmp $outputDir/knownsites/gatk/gatk.HC.vcf
$jumpGatk rm $outputDir/knownsites/gatk/gatk.HC.vcf.idx
$jumpGatk $software->{GATK} SelectVariants -O $outputDir/knownsites/gatk/gatk.SNP.vcf --variant $outputDir/knownsites/gatk/gatk.HC.vcf -select-type SNP
$jumpGatk $software->{GATK} SelectVariants -O $outputDir/knownsites/gatk/gatk.INDEL.vcf --variant $outputDir/knownsites/gatk/gatk.HC.vcf -select-type INDEL

#Keep the intersection of samtools and GATK
$jumpIntersection perl $intersectionScript $outputDir/knownsites/samtools $outputDir/knownsites/gatk $outputDir/knownsites
$jumpIntersection $software->{GATK} IndexFeatureFile -F $outputDir/knownsites/SNP.vcf
$jumpIntersection $software->{GATK} IndexFeatureFile -F $outputDir/knownsites/INDEL.vcf

__COMMAND__
    print $handle $command;
}

sub getSoftwarePath {
    my $file = shift;
    open my $in, '<', $file or die "Error: cannot open file '$file': $!";
    my %software;
    while (<$in>) {
        s/[\r\n]+//;
        next if /\A#/ or ! $_;
        my ($name, $path) = split /\s*=\s*/;
        die "Error: cannot find '$path', please check the file '$file'.\n" if ! -f $path;
        $software{$name} = $path;
    }
    close $in;
    return \%software;
}
