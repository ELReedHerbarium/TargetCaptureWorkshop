###### tags: `target_capture_workshop`

# SESSION 6: Variant Calling


Single Nucleotide Polymorphisms (SNPs), also known as variants are single base pairs of DNA that can be used to detect differences across species and even populations of organisms. Using The Genome Analysis Toolkit (GATK) we will detect or call variants against a reference sample for our set of genes. 


## Variant Calling with Target Capture Data

We will be using a combination of three softwares to call heterozygous variants within a sample: 

* [GATK](https://gatk.broadinstitute.org/hc/en-us) - Genome Analysis Tool Kit
* [bwa](http://bio-bwa.sourceforge.net/bwa.shtml) - Burrows-Wheeler alignment 
* [samtools](http://www.htslib.org/) - toolkit for working with sequence alingment files

### Create a mamba environment

```
mamba create -n variantcall gatk4 bwa samtools
```

Once the installation is complete, activate the mamba environment:

```
mamba activate variantcall
```

### Create A Reference File

Calling variants always requires sequencing reads and a reference genome. Because targeted sequencing frequently uses non-model organisms, there often is not a genome available. We will create a reference file for a single sample using the "supercontig" output of HybPiper (exons and flanking non-coding regions). If you did not run `hybpiper assemble` with the `--run_intronerate` flag in Session 3, go back and do that now before continuing.

#### Action 3.1.1
Concatenate all supercontigs into one single reference file for your species. If the output directory for your HybPiper run is called `prefix` the supercontigs can be recovered with:

```
cat prefix/*/prefix/sequences/intron/*_supercontig.fasta > prefix.supercontigs.fasta
```

*Note: `prefix.supercontigs.fasta` will be used as an input on command line along with* `samplename` below.


### Call SNP Variants Against Reference

The following section will make use of a "Shell Script" to execute many commands on the same set of input files.

Copy the `variantcall.sh` script from: `/home/joh97948/TC_workshop/variantcall.sh`

```

```

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
 
You can use the SFTP/SSH client BitVise or CyberDuck to transfer the following files to your computer: 

- reference supercontigs `prefix.supercontigs.fasta`
- index file for supercontigs `prefix.supercontigs.fasta.fai`
- merged, deduplicated alignment file `prefix.marked.bam`
- index of BAM file `prefix.marked.bam.bai`
- VCF file of filtered SNPs `prefix.snpfiltered.vcf`
- indexed VCF file `prefix.snpfiltered.vcf.idx`


### IGV

The [Integrative Genome Viewer online app](https://igv.org/app/) can be used to view the reference supercontigs and filtered SNPs. You will be able to see the depth of sequences and the locations of variants within the sample (heterozygous sites).

Load the `fasta` file as the "Genome" and the VCF and BAM files as a "Track". Navigate from gene to gene to explore the variation within the sample at different genes.

An IGV viewer guide can be found here: https://software.broadinstitute.org/software/igv/UserGuide
