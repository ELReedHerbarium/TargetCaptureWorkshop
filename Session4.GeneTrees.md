###### tags: `target_capture_workshop`

# SESSION 4: Gene Trees

[TOC]

In this session we will be using the sequences retrieved from Session 3 and running a few commands to be able to arrange the sequences into a gene tree. 

### Download Data

In the previous session, we ran HybPiper to extract data from a single sample. Most projects will have many samples, but in the interest of time we have prepared a set of unaligned FASTA sequences for each of the Angiosperms353 genes, containing 50 taxa from the Kew Tree of Life Explorer.

You can download the data by cloning the github repository for this workshop:

`git clone https://github.com/ELReedHerbarium/TargetCaptureWorkshop.git`

Copy the unaligned sequences to your current directory and uncompress the file:

`cp TargetCaptureWorkshop/data/workshop_paftol_dna.tar.gz .`

`tar -zxf workshop_paftol_dna.tar.gz`

The `workshop_paftol_dna` directory will have 353 FASTA files. Each file will have sequeces from the same 50 taxa. Not all files will have 50 sequences, as we did not recover all genes for all samples.

We will need to make a list of gene names for the next steps. Here you will learn some quick UNIX tricks working with text:

1. Go to the `workshop_paftol_dna` directory and type `ls *.fasta`. This will produce a list of file names like `6457.paftol.fasta`  `7572.paftol.fasta`, etc.
2. We can split the text output of `ls` using the `cut` command and specifying the `.` as a delimeter and choosing the first field:
```
ls *.fasta | cut -f 1 -d '.' 
```
3. We can redirect the output of the command above to a file called `genelist.txt`:

```
ls *.fasta | cut -f 1 -d '.' > genelist.txt
```

Congratulations, you've just written a "pipeline" - a set of UNIX commands that uses the `|` and `>` characters to direct output of one program to another!


### Create mamba environment



We will prepare a `mamba` environment as you did in the last session to install the three softwares we will use today: `mafft`, `trimal`, and `iqtree`

```
mamba create -n phylo mafft trimal iqtree parallel
```

Activate your new mamba environment with `mamba activate phylo`


### MAFFT

When recovering sequences using HybPiper, as we did last session, the sequences are recovered *unaligned*:

![](https://i.imgur.com/tP5dKhs.png)

To properly infer a phylogeny, we wll need to *align* the sequences, so that each column represents homology (shared ancestry). 


The Multiple Alignment using Fast Fourier Transform (MAFFT) pipeline is a command used to create multiple sequence alignments of amino acids or nucleotide sequences from raw sequences. This is an important tool for analysis of the sequences.


First, create a new directory for your MAFFT output:

`mkdir MAFFT`

Here is the base command for running MAFFT

`mafft --preservecase --maxiterate 1000 --localpair input.fasta > output.fasta`

`--preservecase` makes sure that all letters in the sequence stay uppercase

`--maxiterate` is the maximum number of iterative refinement you are allowing

`--localpair` is a more careful iteration that might be slower but is more accurate 

The above command would allow you to align one fasta gene file at a time. We're now working with multiple gene files, so it will help to run our alignments in parallel instead of one at a time. We'll use [GNU parallel](https://www.gnu.org/software/parallel/), which has a peculiar syntax.

It will be much faster to run parallel as follows: 

```
 parallel --eta "mafft --preservecase --maxiterate 1000 --localpair {}.paftol.fasta > MAFFT/{}.mafft.fasta" :::: genelist.txt
```

The way the `parallel` syntax works is that each item in `genelist.txt` will be inserted wherever there is `{}`:

- `{}.workshop.dna.fasta` - the input, unaligned sequence file
- `MAFFT/{}.mafft.fasta` - the output, aligned sequence file

The `--eta` flag for `parallel` gives an estimation of how long until all of the genes finish.

Here is what an aligned sequence file might look like:

![](https://i.imgur.com/ERJfmqc.png)



### Trimal
Aligned sequences are ready for phylogenetic inference, but as you can see above there are multiple parts of the alignment that are not present in every sample. If these *gaps* are too frequent, it can negatively affect the inference of the gene tree.



Trimal is a command that is used to automatically remove any illegitimate or poorly aligned sections from the multiple aligned sequences. Trimal is used to make the alignments an ideal size for placing them on a tree.

First, make a directory for Trimal output:

`mkdir trimal`

Here is the base command for running trimal

`trimal -in <inputfile> -out <outputfile> -(other options)`

You will need your input file (-in) that can be in several different formats. (cluster,fasta,nexus,phylip,NBRF/PIR)

Then you will need to tell the command where you want your output file (-out) to go and what name you want the file to be called

Example

`trimal -in MAFFT/genename.mafft.fasta -out trimal/genename.trimal.fasta -gt .5`

`-gt` is an option that tells trimal how big of a gap is allowed in that fraction of the sequence (in this case, 50%).

To run on all of the genes, try incorporating the parallel command and your genelist. 

Use the space below to show how you will run the trimal command with parallel. 

```
```

A trimmed alignment might look like this:

![](https://i.imgur.com/5IF5GST.png)


### IQ Tree
IQ tree command will take your input of multiple sequence alignment and will reconstruct a phylogeny that is best explained by your input data.

We can now start to reconstruct a maximum-likelihood tree from the alignments you received. Make sure you are in the same directory as your trimal file.

First, change to the trimal output directory:

`cd trimal`

Next, run IQTree on the trimal output file:

`iqtree -s geneName.trimal.fasta`

`-s` is an option to specify the name of the alignment. This command is always required in order for IQtree to work. 

You will receive three output files after running IQtree. 

    geneName.trimal.fasta.iqtree
    
    geneName.trimal.fasta.treefile
    
    geneName.trimal.fasta.log
    
.iqtree is the main file that you should look at to see the calculations. It also contains a view of the final tree in text format.

.treefile is the evolutionary tree that can be viewed in programs like [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) and [iTOL](https://itol.embl.de/).

.log is a file of the entire command run.

Here's what a tree file for your gene tree might look like:
![](https://i.imgur.com/tr1AOtp.png)


#### Determing what model is best
IQtree supports a wide range of models for DNA, protein, codons, binary and nonbinary alignments. If you are not sure which model to use for your file you can run ModelFinder to figure out which is best.

`-m MFP` is the option to specify the model name to use during the analysis. The special MFP key word stands for ModelFinder Plus, which tells IQ-TREE to perform ModelFinder and the remaining analysis using the selected model. 

#### Assessing branch support with bootstrap estimate
IQtree presents an ultrafast bootstrap estimate that is faster than the standard procedure and provides relatively unbiased bootstrap support estimates. 

To run UF boot

`iqtree -s geneName.trimal.fasta -m MFP -B 1000 --redo`

`-B` specifies the number of bootstrap replicates where 1000 is the minimum number recommended.

`--redo` tells IQTree that it's ok to rerun the analysis and overwrite the previous files

Use the space below to show how you will run IQtree in parallel

```
```

A tree with bootstrap labels might look like this:

![](https://i.imgur.com/QDkFkEn.png)


