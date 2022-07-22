#!/bin/bash
set -eo pipefail

## Usage: bash variantcall.sh Sporobolus_cryptandrus_ref.fasta GUMOseq-0188
## Where .fasta is reference sequence + ' ' + samplename 

# Set variables
reference=$1
prefix=$2
read1fn="$prefix".R1.paired.fastq.gz
read2fn="$prefix".R2.paired.fastq.gz

if [ ! -f ../*.dict ]
then gatk CreateSequenceDictionary -R $reference
fi 

bwa index $reference
samtools faidx $reference

#Align read files to reference sequence and map
bwa mem $reference $read1fn $read2fn | samtools view -bS - | samtools sort - -o "$prefix.sorted.bam"
gatk FastqToSam -F1 $read1fn -F2 $read2fn -O $prefix.unmapped.bam -SM $prefix.sorted.bam

#Replace read groups to mapped and unmapped bam files using library prep and sequencing information
gatk AddOrReplaceReadGroups -I  $prefix.sorted.bam -O $prefix.sorted-RG.bam -RGID 2 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM $prefix
gatk AddOrReplaceReadGroups -I  $prefix.unmapped.bam -O $prefix.unmapped-RG.bam -RGID 2 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM $prefix

#Combine mapped and unmapped BAM files
gatk MergeBamAlignment --ALIGNED_BAM $prefix.sorted-RG.bam --UNMAPPED_BAM $prefix.unmapped-RG.bam -O $prefix.merged.bam -R $reference

#Remove duplicate sequences
gatk MarkDuplicates -I $prefix.merged.bam -O $prefix.marked.bam -M $prefix.metrics.txt
samtools index $prefix.marked.bam

#Create VCF 
gatk HaplotypeCaller -I $prefix.marked.bam -O $prefix.vcf -R $reference

#Select only SNPs
gatk SelectVariants -V "$prefix".vcf -R $reference -select-type-to-include SNP -O "$prefix".SNPall.vcf

#Hard filter on variants
gatk VariantFiltration -R $reference -V "$prefix".SNPall.vcf --filter-name "hardfilter" -O "$prefix".snp.filtered.vcf --filter-expression "QD < 5.0 && FS > 60.0 && MQ < 45.0 && MQRankSum < -12.5 && ReadPosRankSum < -8.0"


#index VCF
bgzip $prefix.snp.filtered.vcf
tabix $prefix.snp.filtered.vcf.gz

#Remove intermediate files
rm $prefix.sorted.bam $prefix.unmapped.bam $prefix.merged.bam $prefix.unmapped-RG.bam $prefix.sorted-RG.bam

