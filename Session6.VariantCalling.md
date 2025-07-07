###### tags: `target_capture_workshop`

# SESSION 6: Variant Calling


Single Nucleotide Polymorphisms (SNPs), also known as variants are single base pairs of DNA that can be used to detect differences across species and even populations of organisms. Using The Genome Analysis Toolkit (GATK) we will detect or call variants against a reference sample for our set of genes. 


## Variant Calling with Target Capture Data

We will be using a combination of three softwares to call heterozygous variants within a sample: 

* [GATK](https://gatk.broadinstitute.org/hc/en-us) - Genome Analysis Tool Kit
* [bwa](http://bio-bwa.sourceforge.net/bwa.shtml) - Burrows-Wheeler alignment 
* [samtools](http://www.htslib.org/) - toolkit for working with sequence alingment files
* [whatshap](https://whatshap.readthedocs.io/en/latest/) - for phasing sequences into haplotypes

### Create A Reference File

Calling variants always requires sequencing reads and a reference genome. Because targeted sequencing frequently uses non-model organisms, there often is not a genome available. We will create a reference file for a single sample using the "supercontig" output of HybPiper (exons and flanking non-coding regions). If you did not run `hybpiper assemble` in Session 3, go back and do that now before continuing.


#### Action 3.1.1

If you are in the new `variantcall` conda environment, deactivate it and re-activate the `hybpiper` environment:

```
conda deactivate
```

```
conda activate hybpiper
```

The `hybpiper retrieve_sequences` command can concatenate all supercontigs into one single reference file for your species. If the output directory for your HybPiper run is called `prefix` the supercontigs can be recovered with:

```
hybpiper retrieve_sequences --single_sample_name prefix -t_dna mega353.fasta supercontig
```


This will create a file in your current directory called `prefix.supercontig.fasta`



### Create a conda environment

```
conda create -n variantcall gatk4 bwa samtools whatshap
```

Once the installation is complete, activate the conda environment:

```
conda activate variantcall
```


### Call SNP Variants Against Reference

The following section will make use of a "Shell Script" to execute many commands on the same set of input files.


Copy the `variantcall.sh` script from the GitHub repository: `TargetCaptureWorkshop/scripts/variantcall.sh`



#### Modify shell script

The shell script uses the variable `$prefix` throughout the script to identify the name of the sample. For example, these lines will identify the read files used for mapping reads to the supercontig reference file:
```
read1fn="$prefix"_1P.fastq
read2fn="$prefix"_2P.fastq
```

You may need to modify these lines to match the file names of your trimmed read files. For example, if your files are `ERR0240004.trimmed.R1.fastq` and `ERR0240004.trimmed.R1.fastq` you should modify the two lines above to read:

```
read1fn="$prefix".trimmed.R1.fastq
read2fn="$prefix".trimmed.R2.fastq
```

In this example, `ERR0240004` is the prefix. Use `nano` to modify your shell script. 

#### Run the `variantcall.sh` script

From the command line, you will need to specify the file containing your supercontigs, and also the `$prefix`. Use `bash` to tell the system how to interpret the commands within the `variantcall.sh` script

```
bash variantcall.sh prefix.supercontigs.fasta prefix
```

### What's going on in the script?

The `variantcall.sh` script uses:

- ```bwa mem``` to map paired-end reads to supercontigs.
- Unmapped reads are merged with mapped reads using `gatk` 
- Duplicate reads (identified by having same start and end location) are removed with `gatk RemoveDuplicates`. 
- `gatk HaplotypeCaller` is used to call variants within the sample (heterozygous sites)
- `gatk SelectVaraints` selects only single-nucleotide-polymorphisms (SNPs)
- `gatk VariantFiltration` filters variants based on a hard filtering scheme to remove SNPs not 
- At the end the script removes intermediate BAM files to save space.


#### Action 3.2.2
 
Use the Data browser in the Discovery Environment on Cyverse to transfer the following files to your computer: 

- reference supercontigs `prefix.supercontigs.fasta`
- index file for supercontigs `prefix.supercontigs.fasta.fai`
- merged, deduplicated alignment file `prefix.marked.bam`
- index of BAM file `prefix.marked.bam.bai`
- VCF file of filtered SNPs `prefix.snpfiltered.vcf.gz`
- indexed VCF file `prefix.snpfiltered.vcf.idx`


### IGV

The [Integrative Genome Viewer online app](https://igv.org/app/) can be used to view the reference supercontigs and filtered SNPs. You will be able to see the depth of sequences and the locations of variants within the sample (heterozygous sites). An IGV viewer guide can be found here: https://software.broadinstitute.org/software/igv/UserGuide


Under Genome select "Upload File" and select both the `fasta` file and the `fai` files. 

Next, under Tracks, select Load File and upload the other four files (VCF/index and BAM/index). 

In the top bar where it says `all` you can navigate from gene to gene to explore the variation within the sample at different genes.

![](https://i.imgur.com/IcUdIA8.png)

In the image above, the sequence of one gene is displayed as a "ribbon" with four colors representing A, C, T, and G.

The next "track" is the read depth histogram - some areas of the gene are represented by more reads than others. This is expected for target capture data; we especially expect the depth in regions flanking the targeted exons to drop to zero.

Below the histogram is the read track, where all deviations from the reference are indicated by colors. Finally, at the bottom there are the locations of the filtered SNPs in the VCF track. 

Navigate to a few different genes and answer the following questions:

* Which genes have the most heterozygous sites?
* How does read depth differ from gene to gene?
* Why aren't all of the variants in the read track counted as SNPs in the VCF track?

## Read-Backed Phasing

The VCF created in the previous step treats every variant site (SNP) independently. However, one of the advantages of target capture data is the potential for long haplotype sequences. Inferring which variants are connected across variant sites is known as *phasing* and can be done with WhatsHap, using the variant file (VCF), reads (BAM), and reference sequence (FASTA):

```
whatshap phase -o prefix.phased.vcf prefix.snp.filtered.vcf.gz prefix.marked.bam -r prefix_supercontig.fasta
```

WhatHap also has a nice feature called `haplotag` that will help visualize the reads based on which phase they belong to. First, we need to index the phased VCF:

```
bgzip prefix.phased.vcf
tabix prefix.phased.vcf.gz
```

Next, run `whatshap haplotag`:

```
whatshap haplotag prefix.phased.vcf.gz prefix.marked.bam -r prefix_supercontig.fasta -o prefix.haplotag.bam
```

Now index the new haplotag file:

```
samtools index prefix.haplotag.bam
```

Download the new `haplotag.bam` and index `haplotag.bam.bai` files to your computer and add them to IGV. Next to the new haplotag track in IGV, click on the gear and in the "Color By" selection select "tag". In the dialog that appears, type `HP`

The reads will now be colored based on their haplotype phase. You can now see how SNPs are connected within haplotypes at each gene:

![](https://i.imgur.com/ruxYud1.png)

In the gene above, all positions are connected by reads, so the entire sequence can be phased. This may not be the case for all genes. When the gene cannot be connected, it is separated in to "Phase Blocks" - SNPs within each phase block can be joined into haplotypes, but they cannot be joined across haplotypes.

We can view the phase blocks in IGV by repeating the "Color By... tag" option but using `PS`. Here's what that looks like for another gene:

![](https://i.imgur.com/MjC6NE1.png)

The reads in green contribute to a phase block on the left and the reads in orange contribute to a phase block on the left.

* What are two reasons why we might not be able to phase all SNPs in a gene using target capture data?
