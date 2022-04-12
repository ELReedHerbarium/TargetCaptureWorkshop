###### tags: `target_capture_workshop`

# SESSION 3: HybPiper
[TOC]


## 3.0 Managing Programs with Conda

It can be a challenge to maintain programs on a remote system, especially if the programs are only distributed as source code and need to be compiled on each system separately. In Bioinformatics, often different workflows require different versions of the same software, and these can also be very difficult to maintain.

You have already seen one solution to this issue (Singularity), and this week we will introduce another solution: `conda`. [From the Conda website](https://docs.conda.io/en/latest/):

> Conda is an open source package management system and environment management system that runs on Windows, macOS and Linux. Conda quickly installs, runs and updates packages and their dependencies. Conda easily creates, saves, loads and switches between environments on your local computer. It was created for Python programs, but it can package and distribute software for any language.

### Installing Conda

If you have not worked with Conda on HPCC before, you will need to run the basic installation script. We will be using a version of conda known as "Miniconda" that has a small number of packages included. In a web browser navigate to the [MiniConda installation page](https://docs.conda.io/en/latest/miniconda.html). 

In your Terminal application, login to the TTU HPCC. 

Right-click on the `Miniconda3 Linux 64-bit` link and select "Copy Link Address"

In your remote login terminal use `wget` to download the Miniconda installation script:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Once the script downloads, you will need to run the script:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the prompts on the screen during installation. The default locations for folders will be fine for this tutorial.

When the installation is finished, you will need to logout from HPCC and log back in (because the script modified your `.bashrc` file)

### Setting up Conda

Conda can install packages from a variety of sources - there are default packages and additional packages in *repositories*. When you install packages with `conda` you need to tell `conda` where to look. For bioinformatics packages, Bioconda is a common respository. You need to set up conda to look in each of these repositories. **Run these commands in order**

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Installing via Conda

Installing programs and packages with `conda` is (usually) as simple as `conda install packagename`. For this tutorial we will need to install a second package manager `mamba`:

`conda install mamba`

Follow the prompts on screen. This should take a couple of minutes.

### Creating a Mamba environment

Once you have installed Miniconda, you are ready to create your first environment. Each environment will contain a separate version of Python and any other software you choose to install. 

This week, we will be working with the bioinformatics pipeline `hybpiper` so we will name our environment accordingly:

```
mamba create -n hybpiper
```

Follow the on-screen prompts as the environment is created. Once the evnironment is created, activate it with:

```
mamba activate hybpiper
```

Note that your command prompt has changed with `(hybpiper)` at the beginning to indicate what environment you are in. If you ever need to get out of the envirnoment, use `mamba deactivate`.

### Installing Software in a Conda Environment

Make sure you have the HybPiper environment active by checking that your command prompt begeins with `(hybpiper)` Then install `hybpiper` with this command:

```
mamba install -c chrisjackson-pellicle hybpiper
```

Here, the `-c` is telling `mamba` to look for a package called `hybpiper` within a specific repository (in this case belonging to `chrisjackson-pellicle`, one of the developers of HybPiper). The package repository has a list of required package (prerequisites) that it will also install. 



## 3.1 Running HybPiper
The [HybPiper](https://github.com/mossmatters/HybPiper/blob/master/README.md) bioinformatics pipeline is designed to assemble target capture data, which is when extracted DNA is probed for genes regions of interest. Hybpiper is composed of numerous python scripts that use bioinformatics tools (like BLAST, EXONERATE, and GNU Parallel), to make assembling sequences easier for the user.

HybPiper runs the following tools in this order:
1. Assigns high-throughput reads to target genes using BLASTx or BWA
2. Reads are distributed into directories, and then assembled using SPAdes
3. Creates output FASTA files; one containing CDS portion of the sample for each target region, and the other containing the translated protein sequence.

### Interactive Job With Multiple Processors

HybPiper will assemble sequences from each targeted gene separately. For efficiency, we can make use of multiple processors on the nodes of the HPCC. Previously, when running `interactive` you were only requesting one processor, but to request multiple processors:

```
interactive -p nocona -c 4
```


### Accessing your Target File

First, we will need a FASTA file containing the target genes selected for the dataset. Your target file, if used with a specific kit, can be downloaded from online and FTPed using CyberDuck into your terminal. If you have made your own probes, you will have to create this file yourself.

The reads you downloaded in the previous section were generated with the [Angiosperms353 targeted sequencing kit](https://arborbiosci.com/genomics/targeted-sequencing/mybaits/mybaits-expert/mybaits-expert-angiosperms-353/), so we will use a file 

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

> `hybpiper assemble -t [Target File].fasta -r [Sample File].fastq --prefix [Output filename] --bwa --cpu N --targetfile_ambiguity_codes NMRWSYKVHD`
 
`assemble` is the python code which will run this portion of HybPiper.
 
The target file (`-t`) is the file we have located and moved into your working directory from before.

The sample file fastq file(s) (`-r`) indicates which sample we will be running. You should use the wildcard shortcut we used previously, and have the python script run the forward and reverse reads simultaneously. *Hint: try using the wildcard after the R in your filename.*

 `--bwa` should be used when neucleotides are used in the target files or input files.
 
 `--prefix` will control the name of the output files.
 
`--cpu` specifies the number of processors used by HybPiper. This should match the number of processors you requested in the `interactive` command
 
`--targetfile_ambiguity_codes` lets HybPiper know that this is a target file that contains [DNA ambiguity codes](https://droog.gs.washington.edu/parc/images/iupac.html).

#### Action 3.2
*Try running `hybpiper assemble` by yourself, with the correct filenames.*
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
A supercontig is a ordered set of contigs, creating a portion of the genome. A singular contig is a continuous length of sequence that we are very sure the order of bases is. A supercontig is a larger portion of the genome reconstructed using these smalled contigs, creating a sequence with a few gaps. To recover the supercontigs we will need to re-run `hybpiper assemble` but now with some additional flags:

```
--run_intronerate --no-blast --no-distribute --no-assemble
```

The `--run_intronerate` flag recovers regions flanking the targeted exons, while the various `--no-` commands skip the earlier parts of the workflow.

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
usage: hybpiper stats [-h] [--seq_lengths_filename SEQ_LENGTHS_FILENAME] [--stats_filename STATS_FILENAME] targetfile {dna,DNA,aa,AA} {gene,GENE,supercontig,SUPERCONTIG} namelist
hybpiper stats: error: the following arguments are required: targetfile, targetfile_sequence_type, sequence_type, namelist
```

Therefore we need to specify the targetfile (`mega353.fasta`), the sequence type for the target file (`dna`), the sequence type to summarize (`gene`) and a namelist that contains the 

Use `nano` to create a file called `namelist.txt` that contains the directory name for the sample generated by HybPiper (specified using `--prefix`).

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

Since this command will only make sense if multiple samples are run through HybPiper, you will be given the location of `seq_lengths.tsv` run from multiple samples.


#### Action 3.6

Construct a `hybpiper recovery_heatmap` command for the given file location.

```



```






