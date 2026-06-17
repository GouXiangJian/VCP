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
my $VERSION = 'VCP, step5 v1.3 (2026-06-17)';

#get the home directory of VCP
my $vcp_path = abs_path $0;
$vcp_path =~ s/\\/\//g;
$vcp_path =~ s/\/[^\/]+\z//;
my $softwarePath = "$vcp_path/required_software.txt";
my $filterScript = "$vcp_path/script/variantFilter.pl";

#set default options
my $population = "population_genotype";
my $filterMAF = 0.05;
my $filterMR = 0.2;
my $refGenomeFa = "/public/home/xjgou/genome/genome.fa";
my $outputDir = "output";
my $queue = "smp";
my $thread = 10;
my $memory = 100;
my $version;
my $help;

#get options from command line
GetOptions(
    'population=s'       => \$population,
    'filterMAF|ff=f'     => \$filterMAF,
    'filterMR|fr=f'      => \$filterMR,
    'refGenomeFa=s'      => \$refGenomeFa,
    'outputDir=s'        => \$outputDir,
    'queue=s'            => \$queue,
    'thread=i'           => \$thread,
    'memory=i'           => \$memory,
    'version+'           => \$version,
    'help+'              => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
####################################################################################################
Name: VCP - Variant Calling Pipeline (step5)

Function: variant filtering, including: GatherVcfs, VariantFiltration, SelectVariants, and Tassel5

Usage: perl step5.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  #Options for filtering:
  -p  | -population   <STRING> : the population name. (default: population_genotype)
  -ff | -filterMAF    <FLOAT>  : the minor allele frequency threshold. (default: 0.05)
  -fr | -filterMR     <FLOAT>  : the missing rate threshold. (default: 0.2)

  #Options for path:
  -r  | -refGenomeFa  <STRING> : the reference genomic fasta file. (default: /public/home/xjgou/genome/genome.fa)
  -o  | -outputDir    <STRING> : set a directory for storing output information. (default: output)

  #Options for resources:
  -q  | -queue        <STRING> : set the queue to use. (default: smp)
  -t  | -thread       <INT>    : set the number of threads to use. (default: 10)
  -m  | -memory       <INT>    : set the size of memory to use. (default: 100 [100GB])

  #Options for other:
  -v  | -version               : show the version information.
  -h  | -help                  : show the help information.
####################################################################################################

__GUIDE__

#output version and help information
die "$VERSION\n" if $version;
die $usage if $help;

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
my $input = join ' ', map {"-I $outputDir/GenotypeGVCFsOut/$_.vcf"} @names;

#generate a scheduling script for filtering, meanwhile, generate a comprehensive scheduling script
my $lsfDir = "step5.lsf";
system "mkdir -p $lsfDir";
open my $oTOTAL, '>', 'step5.sh';
print $oTOTAL "bsub < $lsfDir/1.lsf\n";
close $oTOTAL;
open my $oEACH, '>', "$lsfDir/1.lsf";
lsfInfo($oEACH, $queue, $thread, $memory, $input, $outputDir, $softwarePath, $population, $filterScript, $filterMAF, $filterMR);
close $oEACH;

#create a subroutine to write all command into lsf script
sub lsfInfo {
    my ($handle, $queue, $thread, $memory, $input, $outputDir, $softwarePath, $population, $filterScript, $filterMAF, $filterMR) = @_;
    my $software = getSoftwarePath($softwarePath);
    my $command = <<__COMMAND__;
#!/bin/bash

#BSUB -q $queue
#BSUB -n $thread
#BSUB -J step5
#BSUB -o step5.out
#BSUB -e step5.err

#run GatherVcfs
$software->{GATK} GatherVcfs $input -O $outputDir/GenotypeGVCFsOut/$population.vcf

#run VariantFiltration
mkdir -p $outputDir/final
$software->{GATK} VariantFiltration -V $outputDir/GenotypeGVCFsOut/$population.vcf -O $outputDir/final/$population.HC.vcf --cluster 4 --window 10 --mask-extension 3 --filter-name "lowMQ" --filter "MQ < 40.0" --filter-name "lowDP" --filter "DP < 8.0" --filter-name "LowQual" --filter "QUAL < 20" --filter-name "lowQD" --filter "QD < 2.0" --filter-name "lowReadPosRankSum" --filter "ReadPosRankSum < -8.0" --filter-name "highFS" --filter "FS > 60.0" --filter-name "lowMQRankSum" --filter "MQRankSum < -12.5"

#run bcftools
$software->{BCFTOOLS} view -f PASS $outputDir/final/$population.HC.vcf > $outputDir/final/tmp
mv $outputDir/final/tmp $outputDir/final/$population.HC.vcf
rm $outputDir/final/$population.HC.vcf.idx

#run SelectVariants
$software->{GATK} SelectVariants -O $outputDir/final/$population.SNP.vcf --variant $outputDir/final/$population.HC.vcf -select-type SNP
$software->{GATK} SelectVariants -O $outputDir/final/$population.INDEL.vcf --variant $outputDir/final/$population.HC.vcf -select-type INDEL

#run SortGenotypeFilePlugin
perl $software->{TASSEL5} -Xmx${memory}g -SortGenotypeFilePlugin -inputFile $outputDir/final/$population.SNP.vcf -outputFile $outputDir/final/$population.SNP.sort.vcf -fileType VCF
perl $software->{TASSEL5} -Xmx${memory}g -SortGenotypeFilePlugin -inputFile $outputDir/final/$population.INDEL.vcf -outputFile $outputDir/final/$population.INDEL.sort.vcf -fileType VCF

#get final hmphap file
perl $software->{TASSEL5} -Xmx${memory}g -fork1 -vcf $outputDir/final/$population.SNP.sort.vcf -export $outputDir/final/$population.SNP -exportType Hapmap -runfork1
perl $software->{TASSEL5} -Xmx${memory}g -fork1 -vcf $outputDir/final/$population.INDEL.sort.vcf -export $outputDir/final/$population.INDEL -exportType Hapmap -runfork1

#filter based on MAF and MR (only keep bi-allelic variants)
perl $filterScript $outputDir/final/$population.SNP.hmp.txt $filterMAF $filterMR $outputDir/final/$population.SNP
perl $filterScript $outputDir/final/$population.INDEL.hmp.txt $filterMAF $filterMR $outputDir/final/$population.INDEL

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
