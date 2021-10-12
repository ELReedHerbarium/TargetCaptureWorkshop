# SESSION 3: HybPiper
[TOC]

## 3.1 Running HybPiper
The [HybPiper](https://github.com/mossmatters/HybPiper/blob/master/README.md) bioinformatics pipeline is designed to assemble target capture data, which is when extracted DNA is probed for genes regions of interest. Hybpiper is composed of numerous python scripts that use bioinformatics tools (like BLAST, EXONERATE, and GNU Parallel), to make assembling sequences easier for the user.

HybPiper runs the following tools in this order:
1. Assigns high-throughput reads to target genes using BLASTx or BWA
2. Reads are distributed into directories, and then assembled using SPAdes
3. Creates output FASTA files; one containing CDS portion of the sample for each target region, and the other containing the translated protein sequence.

### Accessing your Target File

To get started, we need to either locate and move files, or create them. First, we will need a fasta file containing the target genes selected for the dataset. For example, when using the Angiosperms353 Target capture kit, we should use the target file called Mega353.fasta. Your target file, if used with a specific kit, can be downloaded from online and FTPed using CyberDuck into your terminal. If you have made your own probes, you will have to create this file yourself. Use the space below to show how you will move the target file into your working directory.

```

```

### Creating a Namelist

Second, we will need to create a namelist text file to help the pipeline run smoothly. While this is not required to run the program, it makes running multiple samples easier. Use the `nano` command to create a .txt file that has sample names, one per line.
```

```

### Running HybPiper

With these two files, we are able to run the first HybPiper command, `reads_first.py`. The arguments for this command are as follows:

> `../reads_first.py -b [Target File].fasta -r [Sample File].fastq --prefix [Output filename] --bwa`
> reads_first.py is the python code which will run this portion of HybPiper.
> 
> The target file is the file we have located and moved into your working directory from before.
> 
> The sample file fastq file indicates which sample we will be running. You should use the wildcard shortcut we used previously, and have the python script run the forward and reverse reads simultaneously. *Hint: try using the wildcard after the R in your filename.*
> 
> the output file name should be the name of your sample used in the textfile of all sample names.
> 
> -bwa should be used when neucleotides are used in the target files or input files.
> the -b, -r, and --prefix are flags that organize and direct the command.

Try running `reads_first.py` by yourself, with the correct filenames.
```

```

### Using a While Loop

Hybpiper runs each set of reads independently, meaning that we could continue to run `reads_first.py` for each set of data. However, if you wanted to run an entire plate of samples, this process would become quite time consuming. Lets write a simple script to run all of our samples at once.

`
while read name; 
do ../reads_first.py -b test_targets.fa -r $name*.fastq --prefix $name --bwa
done < namelist.txt
`


A while loop is statement that allows for a command to be repeated. We are able to create this while loop because the namelist.txt file acts as a variable for filenames, in addition to setting the  --prefix flag. While loops and the namelist text file will be utilized later in this tutorial when we want to analyze the HybPiper output data.

### Supercontigs
A supercontig is a ordered set of contigs, creating a portion of the genome. A singular contig is a continuous length of sequence that we are very sure the order of bases is. A supercontig is a larger portion of the genome reconstructed using these smalled contigs, creating a sequence with a few gaps. Hybpiper assembled these supercontigs when we run `reads_first.py` .

### Getting Sequence Lengths

To get a quick visual summary of our data thus far, we can use `get_seq_lengths.py`.

`python ../get_seq_lengths.py [Hybpiper output files].fasta [namelist].txt dna > [Seq Lengths output filename].txt`

Replace the example arguments with your own, and run the command.

```

```
The first line of test_seq_lengths.txt has the names of each gene. The second line has the length of the target gene, averaged over each "source" for that gene. The rest of the lines are the length of the sequence recovered by HybPiper for each gene. If there was no sequence for a gene, a 0 is entered.

We can view the output file by opening it with `less`. The first line of the file has all the target gene names, and the second line has the length of the target gene. The lines following will have each sample per line, showing the length of the sequence recovered for that particular gene. If there were no reads mapped to a particular gene, you will see a 0 in that sample line. If you are familliar with a heatmap, this output file it what you would use to create one in RStudio. Visualizing this data is an important step to determining that quality of your data.

### Getting HybPiper Stats

Getting a summary of statistics for our assembled reads is easy using the `hybpiper_stats.py` 

> This script will summarize target enrichment and gene recovery efficiency for a set of samples. The output is a text file with one sample per line and the following statistics:
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

This is the usage for `hybpiper.py`

`python hybpiper_stats.py [YourSeqLengths].txt [namelist].txt > [Your summary Statistics].txt`

```

```

### Retrieving Sequences

The `retrieve_sequences.py` script will show us which sequeces from which sample mapped to a given gene. This script will take all runs from HybPiper (`reads_first.py`), and retrieves the gene names from the target file, creating a output fasta file for each gene.

python ../retrieve_sequences.py [target files].fasta . [sequence type]

This script can use protein sequences, or animo acid sequences to create the fasta. In the place of [sequence type], you can choose `aa` or `dna`.

```

```
### Cleaning Up

HybPiper creates a lot of files we dont want to clog up the cluster. This step will ensure that we save space. 

python ../cleanup.py [HybPiper Output file]

```

```
## 3.2 Find the Helper script in HybPiper

If you are ever having tr







