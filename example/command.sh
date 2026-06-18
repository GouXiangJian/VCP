#==========The following is the command line to run the example data==========#

#Note: it takes about 30 minutes.

#1. change the current working directory (it needs to be modified!)
vcp_path=/public/home/xjgou/tools/VCP-main
cd $vcp_path/example/genome

#2. build genomic library and create index
genome_prefix=Saccharomyces_cerevisiae.R64-1-1.dna.toplevel
gzip -d $genome_prefix.fa.gz
bwa index -p $genome_prefix $genome_prefix.fa
samtools faidx $genome_prefix.fa
gatk CreateSequenceDictionary -R $genome_prefix.fa -O $genome_prefix.dict

#3. variant calling
cd $vcp_path/example
perl $vcp_path/step1.pl -g $vcp_path/example/raw_reads -rl $vcp_path/example/genome/$genome_prefix -rf $vcp_path/example/genome/$genome_prefix.fa -q normal -t 4 -m 20 && sh step1.sh
perl $vcp_path/step2.pl -g $vcp_path/example/raw_reads -r $vcp_path/example/genome/$genome_prefix.fa -q normal -t 4 -m 20 && sh step2.sh
perl $vcp_path/step3.pl -g $vcp_path/example/raw_reads -r $vcp_path/example/genome/$genome_prefix.fa -q normal -t 4 -m 20 && sh step3.sh
perl $vcp_path/step4.pl -r $vcp_path/example/genome/$genome_prefix.fa -q normal -t 4 -m 20 && sh step4.sh
perl $vcp_path/step5.pl -r $vcp_path/example/genome/$genome_prefix.fa -q normal -t 4 -m 20 && sh step5.sh

#Note: the final genotype files (VCF/HapMap) are saved in 'output/final' directory. If VCP runs normally, 11 files
#      will be generated in 'output/final' directory: 5 SNP files, 5 INDEL files, and 1 SNP+INDEL (named HC) file.
