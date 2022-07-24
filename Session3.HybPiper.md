###### tags: `target_capture_workshop`

# SESSION 3: HybPiper
[TOC]


## 3.0 Creating a Mamba environment

We will be working with the bioinformatics pipeline `hybpiper` so we will create a new environment named `hybpiper`:

```
mamba create -n hybpiper
```

Follow the on-screen prompts as the environment is created. Once the evnironment is created, activate it with:

```
mamba activate hybpiper
OR
conda activate hybpiper
```

Note that your command prompt has changed with `(hybpiper)` at the beginning to indicate what environment you are in. If you ever need to get out of the environment, use `mamba deactivate`.

### Installing Software in a Conda Environment

Make sure you have the HybPiper environment active by checking that your command prompt begins with `(hybpiper)` Then install `hybpiper` and `time` with these commands:

```
conda config --add channels bioconda

mamba install -c chrisjackson-pellicle hybpiper

mamba install time
```

The installation will take a few minutes. Here, the `-c` is telling `mamba` to look for a package called `hybpiper` within a specific repository (in this case belonging to `chrisjackson-pellicle`, one of the developers of HybPiper). The package repository has a list of required package (prerequisites) that it will also install. 



## 3.1 Running HybPiper
The [HybPiper](https://github.com/mossmatters/HybPiper/blob/master/README.md) bioinformatics pipeline is designed to assemble target capture data, which is when extracted DNA is probed for gene regions of interest. Hybpiper is composed of numerous python scripts that use bioinformatics tools (like BLAST, EXONERATE, and GNU Parallel), to make assembling sequences easier for the user.

HybPiper runs the following tools in this order:
1. Assigns high-throughput reads to target genes using BLASTx or BWA
2. Reads are distributed into directories, and then assembled using SPAdes
3. Creates output FASTA files; one containing CDS portion of the sample for each target region, and the other containing the translated protein sequence.


### Accessing your Target File

First, we will need a FASTA file containing the target genes selected for the dataset. Your target file, if used with a specific kit, can be downloaded from online and FTPed using CyberDuck into your terminal. If you have made your own probes, you will have to create this file yourself.

The reads you downloaded in the previous section were generated with the [Angiosperms353 targeted sequencing kit](https://arborbiosci.com/genomics/targeted-sequencing/mybaits/mybaits-expert/mybaits-expert-angiosperms-353/), so we will use the following file 

Use `wget` to download the `mega353.fasta` file from Github: https://github.com/chrisjackson-pellicle/NewTargets/raw/master/mega353.fasta.zip

The file is compressed by `zip` so you will need to `unzip` the file:

```
unzip mega353.fasta.zip
```

#### Action 3.1
 *Use the space below to show how you will move the target file into your working directory.*
```

```


### Running HybPiper

With these two files, we are able to run the first HybPiper command, `hybpiper assemble`. The arguments for this command are as follows:

> `hybpiper assemble -t_dna [Target File].fasta -r [Sample File].fastq --prefix [Output filename] --bwa --cpu N`
 
`assemble` is the python code which will run this portion of HybPiper.
 
The target file (`-t`) is the file you previously downloaded into your working directory, `mega353.fasta`. Our target file is dna code, so we should use the prefix `-t_dna`

The sample file fastq file(s) (`-r`) indicates which sample we will be running. You should use the wildcard shortcut we used previously, and have the python script run the forward and reverse trimmed reads simultaneously. *Hint: try using the wildcard after the R in your filename.*

 `--bwa` should be used when nucleotides are used in the target files or input files.
 
 `--prefix` will control the name of the output files.
 
`--cpu` specifies the number of processors used by HybPiper. For today, you can use 4.


#### Action 3.2
*Try running `hybpiper assemble` by yourself, with the correct filenames (trimmed reads).*
```

```

### HybPiper Output

The output of HybPiper is a directory named with the `--prefix` flag. This directory contains a set of standard subdirectories, one per gene. Within each gene, there are other directories organizing outputs from the different stages of HybPiper. 

For example, to find the sequences recovered for gene `4471` for a sample with the directory name `SampleName`, you would look in: 

`SampleName/4471/SampleName/sequences/FNA/4471.FNA`

#### Action 3.3

Inside the output directory, the file `genes_with_seqs.txt` contains a list of all the genes with recovered sequences. How many genes were recovered for your sample?

Use CyberDuck to access the HybPiper output directory and navigate the directories within each gene. Does your sample have a sequence recovered for the following genes:
* `4471`
* `6969`
* `7628`


### Supercontigs
A supercontig is an ordered set of contigs, creating a portion of the genome. A singular contig is a continuous length of sequence that we are very sure the order of bases is. A supercontig is a larger portion of the genome reconstructed using these smalled contigs, creating a sequence with a few gaps. To recover the supercontigs we will need to re-run `hybpiper assemble` but now with some additional flags:

```
--run_intronerate --start_from exonerate_contigs
```

The `--run_intronerate` flag recovers regions flanking the targeted exons, while the `--start_from` command skips the earlier parts of the workflow.

#### Action 3.4

The output of this command will add a directory `intron` to each gene recovered. Inside this directory is a file with the suffix `_supercontig.fasta`

Identify one of the genes with a recovered sequence and download both the `FNA` file and `supercontig.fasta` file. For example, you would find sequences for gene `4471` in:

`SampleName/4471/SampleName/sequences/FNA/4471.FNA` 

and

`SampleName/4471/SampleName/sequences/intron/4471_supercontig.fasta`

Open both in a text editor: **how different are the lengths of the sequences?**



### Getting HybPiper Stats

Getting a summary of statistics for our assembled reads is easy using the `hybpiper stats`:

```
usage: hybpiper stats [-h]
                      (--targetfile_dna TARGETFILE_DNA | --targetfile_aa TARGETFILE_AA)
                      [--seq_lengths_filename SEQ_LENGTHS_FILENAME]
                      [--stats_filename STATS_FILENAME]
                      {gene,GENE,supercontig,SUPERCONTIG} namelist
```

Therefore we need to specify the targetfile (`-t_dna mega353.fasta`),  the sequence type to summarize (`gene`), and a namelist.

Since we only have one sample, we need to use `nano` to create a file called `namelist.txt`. This file should contains the directory name for the sample generated by HybPiper (specified using `--prefix`).

#### Action 3.5
*Replace the example arguments with your own, and run the `hybpiper stats` command.*
```

```

`hybpiper stats` will generate two files, `seq_lengths.tsv` and `hybpiper_stats.tsv`

The first line of `seq_lengths.tsv` has the names of each gene. The second line has the length of the target gene, averaged over each "source" for that gene. The rest of the lines are the length of the sequence recovered by HybPiper for each gene. If there was no sequence for a gene, a 0 is entered.

We can view the output file by opening it with `less`. The first line of the file has all the target gene names, and the second line has the length of the target gene. The lines following will have each sample per line, showing the length of the sequence recovered for that particular gene. If there were no reads mapped to a particular gene, you will see a 0 in that sample line. 

**How many genes are missing sequences?**


The `hybpiper_stats.tsv` contains a number of useful statistics about HybPiper output.


This script will summarize target enrichment and gene recovery efficiency for a set of samples. The output is a text file with one sample per line and the following statistics:
> 
> * Number of reads
> * Number of reads on target
> * Percent reads on target
> * Number of genes with reads
> * Number of genes with contigs
> * Number of genes with sequences
> * Number of genes with sequences > 25% of the target length
> * Number of genes with sequences > 50% of the target length
> * Number of genes with sequences > 75% of the target length
> * Number of genes with sequences > 150% of the target length
> * Number of genes with paralog warnings

### Making a Heatmap

A heatmap can be a good visualization of target recovery with multiple samples; each column is a gene, each row is a sample, and the shading in each cell represents the percentage of the target gene length recovered:

![](https://github.com/mossmatters/HybPiper/wiki/img/test_dataset_recovery_heatmap.png)

You will only have one sample, but you can still use the `seq_lengths.tsv` run from `hybpiper stats`.


#### Action 3.6

Construct a `hybpiper recovery_heatmap` command for your `seq_lengths.tsv` file.

```



```






