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
my $VERSION = 'VCP, step4 v1.3 (2026-06-17)';

#get the home directory of VCP
my $vcp_path = abs_path $0;
$vcp_path =~ s/\\/\//g;
$vcp_path =~ s/\/[^\/]+\z//;
my $softwarePath = "$vcp_path/required_software.txt";

#set default options
my $refGenomeFa = "/public/home/xjgou/genome/genome.fa";
my $outputDir = "output";
my $jumpGenomicsDBImport = '';
my $jumpGenotypeGVCFs = '';
my $queue = "smp";
my $thread = 10;
my $memory = 100;
my $version;
my $help;

#get options from command line
GetOptions(
    'refGenomeFa=s'            => \$refGenomeFa,
    'outputDir=s'              => \$outputDir,
    'jumpGenomicsDBImport|jm+' => \$jumpGenomicsDBImport,
    'jumpGenotypeGVCFs|jt+'    => \$jumpGenotypeGVCFs,
    'queue=s'                  => \$queue,
    'thread=i'                 => \$thread,
    'memory=i'                 => \$memory,
    'version+'                 => \$version,
    'help+'                    => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Name: VCP - Variant Calling Pipeline (step4)

Function: joint genotyping, including: GenomicsDBImport and GenotypeGVCFs

Usage: perl step4.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for path:
  -r  | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/genome/genome.fa)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)

  #Options for step:
  -jm | -jumpGenomicsDBImport  : no execute GenomicsDBImport
  -jt | -jumpGenotypeGVCFs     : no execute GenotypeGVCFs

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
$jumpGenomicsDBImport = '#' if $jumpGenomicsDBImport;
$jumpGenotypeGVCFs = '#' if $jumpGenotypeGVCFs;

#get the all chr name
open my $iFASTA, '<', $refGenomeFa or die "Error: cannot open file '$refGenomeFa': $!";
my @names;
while (<$iFASTA>) {
    if (/\A>/) {
        my ($name) = /\A>(\S+)/;
        push @names, $name;
    }
}
close $iFASTA;

#generate a scheduling script for each chromosome, meanwhile, generate a comprehensive scheduling script
open my $oTOTAL, '>', 'step4.sh';
my $lsfDir = "step4.lsf";
system "mkdir -p $lsfDir";
foreach my $chr (@names) {
    open my $oEACH, '>', "$lsfDir/$chr.lsf";
    lsfInfo($oEACH, $chr, $queue, $thread, $memory, $outputDir, $softwarePath, $refGenomeFa, $jumpGenomicsDBImport, $jumpGenotypeGVCFs);
    close $oEACH;
    print $oTOTAL "bsub < $lsfDir/$chr.lsf\n";
}
close $oTOTAL;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $chr, $queue, $thread, $memory, $outputDir, $softwarePath, $refGenomeFa, $jumpGenomicsDBImport, $jumpGenotypeGVCFs) = @_;
    my $input = join ' ', map {"-V $_"} grep {/vcf\z/} glob "$outputDir/HaplotypeCallerOut/*";
    my $software = getSoftwarePath($softwarePath);
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -J step4.$chr
#BSUB -o step4.$chr.out
#BSUB -e step4.$chr.err

#run GenomicsDBImport
$jumpGenomicsDBImport mkdir -p $outputDir/GenomicsDBImportOut
$jumpGenomicsDBImport $software->{GATK} --java-options "-Xmx${memory}g -Xms${memory}g" GenomicsDBImport --genomicsdb-workspace-path $outputDir/GenomicsDBImportOut/$chr --batch-size 50 -R $refGenomeFa -L $chr $input

#run GenotypeGVCFs
$jumpGenotypeGVCFs mkdir -p $outputDir/GenotypeGVCFsOut
$jumpGenotypeGVCFs $software->{GATK} --java-options "-Xmx${memory}g -Xms${memory}g" GenotypeGVCFs -R $refGenomeFa -V gendb://$outputDir/GenomicsDBImportOut/$chr -O $outputDir/GenotypeGVCFsOut/$chr.vcf -new-qual -G StandardAnnotation --use-new-qual-calculator

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
