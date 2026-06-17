#!/usr/bin/env perl

#date    : 2026-06-17
#writer1 : Xiangjian Gou (xjgou@mail.hzau.edu.cn)
#writer2 : Haoran Shi (sicau3339@outlook.com)

#load modules
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;

#record version information
my $VERSION = 'VCP, step1 v1.3 (2026-06-17)';

#get the home directory of VCP
my $vcp_path = abs_path $0;
$vcp_path =~ s/\\/\//g;
$vcp_path =~ s/\/[^\/]+\z//;
my $softwarePath = "$vcp_path/required_software.txt";

#set default options
my $germplasmDir = "/public/home/xjgou/raw_reads";
my $refGenomeLib = "/public/home/xjgou/genome/genome";
my $refGenomeFa = "/public/home/xjgou/genome/genome.fa";
my $outputDir = "output";
my $mbq = 20;
my $jumpNGSQC = '';
my $jumpBWA = '';
my $jumpSortSam = '';
my $jumpMarkDuplicates = '';
my $jumpHaplotypeCaller = '';
my $queue = "normal";
my $thread = 6;
my $memory = 30;
my $version;
my $help;

#get options from command line
GetOptions(
    'germplasmDir=s'          => \$germplasmDir,
    'refGenomeLib|rl=s'       => \$refGenomeLib,
    'refGenomeFa|rf=s'        => \$refGenomeFa,
    'outputDir=s'             => \$outputDir,
    'mbq|b=i'                 => \$mbq,
    'jumpNGSQC|jn+'           => \$jumpNGSQC,
    'jumpBWA|jb+'             => \$jumpBWA,
    'jumpSortSam|js+'         => \$jumpSortSam,
    'jumpMarkDuplicates|jm+'  => \$jumpMarkDuplicates,
    'jumpHaplotypeCaller|jh+' => \$jumpHaplotypeCaller,
    'queue=s'                 => \$queue,
    'thread=i'                => \$thread,
    'memory|m=i'              => \$memory,
    'version+'                => \$version,
    'help+'                   => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Name: VCP - Variant Calling Pipeline (step1)

Function: read preprocessing, including: NGSQC, BWA, SortSam, MarkDuplicates, and HaplotypeCaller (non-BQSR)

Usage: perl step1.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for path:
  -g  | -germplasmDir <STRING> : the directory of storing all germplasms. (default: /public/home/xjgou/raw_reads)
  -rl | -refGenomeLib <STRING> : the reference genomic library. (default: /public/home/xjgou/genome/genome)
  -rf | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/genome/genome.fa)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)
  -b  | -mbq          <INT>    : minimum base quality required to consider a base for calling (default: 20)

  #Options for step:
  -jn | -jumpNGSQC             : no execute NGSQC
  -jb | -jumpBWA               : no execute BWA
  -js | -jumpSortSam           : no execute SortSam
  -jm | -jumpMarkDuplicates    : no execute MarkDuplicates
  -jh | -jumpHaplotypeCaller   : no execute HaplotypeCaller (non-BQSR)

  #Options for resources:
  -q  | -queue  <STRING>       : set the queue to use. (default: normal)
  -t  | -thread    <INT>       : set the number of threads to use. (default: 6)
  -m  | -memory    <INT>       : set the size of memory to use. (default: 30 [30GB])

  #Options for other:
  -v  | -version               : show the version information.
  -h  | -help                  : show the help information.

Notes:
  (1) -g specifies the input directory, which should contain multiple subdirectories, each
      representing one sample. For N samples, the directory should contain N subdirectories.
      Each subdirectory must contain two files: dirName-1.fq.gz and dirName-2.fq.gz.

  (2) The bundled version of the required software is recommended, and the Perl module
      'String::Approx' must be installed before using NGSQC.

  (3) All files and directories provided as parameters should be absolute paths.

  (4) Remember to build a genomic library and create index (.fai/.dict) before using the script.
        \$ bwa index -p genome genome.fa
        \$ samtools faidx genome.fa
        \$ gatk CreateSequenceDictionary -R genome.fa -O genome.dict
####################################################################################################

__GUIDE__

#output version and help information
die "$VERSION\n" if $version;
die $usage if $help;

#set whether to execute each step
$jumpNGSQC = '#' if $jumpNGSQC;
$jumpBWA = '#' if $jumpBWA;
$jumpSortSam = '#' if $jumpSortSam;
$jumpMarkDuplicates = '#' if $jumpMarkDuplicates;
$jumpHaplotypeCaller = '#' if $jumpHaplotypeCaller;

#generate a scheduling script for each germplasm, meanwhile, generate a comprehensive scheduling script
open my $oTOTAL, '>', 'step1.sh';
my $lsfDir = "step1.lsf";
system "mkdir -p $lsfDir";
foreach my $dir (glob "$germplasmDir/*") {
    my ($germplasm) = $dir =~ /([^\/]+)\z/;
    open my $oEACH, '>', "$lsfDir/$germplasm.lsf";
    lsfInfo($oEACH, $germplasm, $germplasmDir, $outputDir, $softwarePath, $queue, $thread, $memory, $refGenomeLib, $refGenomeFa, $jumpNGSQC, $jumpBWA, $jumpSortSam, $jumpMarkDuplicates, $jumpHaplotypeCaller, $mbq);
    close $oEACH;
    print $oTOTAL "bsub < $lsfDir/$germplasm.lsf\n";
}
close $oTOTAL;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $germplasm, $germplasmDir, $outputDir, $softwarePath, $queue, $thread, $memory, $refGenomeLib, $refGenomeFa, $jumpNGSQC, $jumpBWA, $jumpSortSam, $jumpMarkDuplicates, $jumpHaplotypeCaller, $mbq) = @_;
    my $software = getSoftwarePath($softwarePath);
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -J step1.$germplasm
#BSUB -o step1.$germplasm.out
#BSUB -e step1.$germplasm.err

#create output directory
mkdir -p $outputDir

#run NGSQC
$jumpNGSQC mkdir -p $outputDir/NgsqcOut
$jumpNGSQC perl $software->{NGSQC} -pe $germplasmDir/$germplasm/$germplasm-1.fq.gz $germplasmDir/$germplasm/$germplasm-2.fq.gz N A -l 70 -s 20 -o $outputDir/NgsqcOut/$germplasm

#run BWA
$jumpBWA mkdir -p $outputDir/BwaOut
$jumpBWA $software->{BWA} mem -t $thread -Y -a -M -R "\@RG\\tID:$germplasm\\tLB:$germplasm\\tSM:$germplasm\\tPL:ILLUMINA" $refGenomeLib $outputDir/NgsqcOut/$germplasm/$germplasm-1.fq.gz_filtered $outputDir/NgsqcOut/$germplasm/$germplasm-2.fq.gz_filtered 2> $outputDir/BwaOut/$germplasm.log | $software->{SAMTOOLS} view -bS - > $outputDir/BwaOut/$germplasm.bam

#run SortSam
$jumpSortSam mkdir -p $outputDir/SortSamOut
$jumpSortSam $software->{GATK} --java-options "-Dsamjdk.compression_level=5 -Xms4000m" SortSam --INPUT $outputDir/BwaOut/$germplasm.bam --OUTPUT $outputDir/SortSamOut/$germplasm.sort.bam --SORT_ORDER coordinate --CREATE_INDEX false --CREATE_MD5_FILE false --VALIDATION_STRINGENCY SILENT

#run MarkDuplicates
$jumpMarkDuplicates mkdir -p $outputDir/MarkDuplicatesOut
$jumpMarkDuplicates $software->{GATK} --java-options "-Dsamjdk.compression_level=5 -Xms4000m" MarkDuplicates --INPUT $outputDir/SortSamOut/$germplasm.sort.bam --OUTPUT $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam --METRICS_FILE $outputDir/MarkDuplicatesOut/$germplasm.metrics --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER coordinate
$jumpMarkDuplicates $software->{SAMTOOLS} index $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam

#run HaplotypeCaller (non-BQSR)
$jumpHaplotypeCaller mkdir -p $outputDir/HaplotypeCallerOut_nonBQSR
$jumpHaplotypeCaller $software->{GATK} --java-options "-Xmx${memory}G" HaplotypeCaller -R $refGenomeFa -I $outputDir/MarkDuplicatesOut/$germplasm.sort.mark.bam -O $outputDir/HaplotypeCallerOut_nonBQSR/$germplasm.g.vcf -ERC GVCF -stand-call-conf 30 -mbq $mbq --native-pair-hmm-threads 20

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
