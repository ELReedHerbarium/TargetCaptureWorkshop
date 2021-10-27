###### tags: `target_capture_workshop`

# SESSION 4: Gene Trees

SESSION 4: Gene Tree
[TOC]

In this session we will be using the sequences retrieved from Session 3 and running a few commands to be able to arrange the sequences into a gene tree. 

### MAFFT
The Multiple Alignment using Fast Fourier Transform (MAFFT) pipeline is a command used to create multiple sequence alignments of amino acids or nucleotide sequences from raw sequences. This is an important tool for analysis of the sequences.

To get started, you will need to go to the directory that has your sequences you previously ran in session 3.

Here is the base command for running MAFFT

`mafft --preservecase --maxiterate 1000 --localpair input.fasta > output.fasta`

`--preservecase` makes sure that all letters in the sequence stay uppercase

`--maxiterate` is the maximum number of iterative refinement you are allowing

`--localpair` is a more careful iteration that might be slower but is more accurate 

Before running this command, it is in your best interest to make a directory named MAFFT so the command will send the data to that directory and your current directory remains less cluttered.

`"mafft --preservecase --maxiterate 1000 --localpair [YourSequenceDirectory]/[yourfilename]_supercontig.fasta > MAFFT/[yourfilename].mafft.fasta" :::: [namelist].txt`



`MAFFT/{}.mafft.fasta` is my new file with the aligned sequences.

### Trimal
Trimal is a command that is used to automatically remove any illegitimate or poorly aligned sections from the multiple aligned sequences. Trimal is used to make the alignments an ideal size for placing them on a tree.

Here is the base command for running trimal

`trimal -in <inputfile> -out <outputfile> -(other options)`

You will need your input file (-in) that can be in several different formats. (cluster,fasta,nexus,phylip,NBRF/PIR)

Then you will need to tell the command where you want your output file (-out) to go and what name you want the file to be called

Example

`trimal -in [DirectorywhereMAFFTfileIs]/$name.mafft.fasta -out [DirectoryWhereyouWantFiletoGo]/$name.trimal.fasta -gt .5 :::: genelist.txt`

`-gt` is an option that tells trimal how big of a gap is allowed in that fraction of the sequence.

Before running this command, it is in your best interest to make a directory named Trimal so the command will send the data to that directory and your current directory remains less cluttered.

Use the space below to show how you will run the trimal command 

```
```


### IQ Tree
IQ tree command will take your input of multiple sequence alignment and will reconstruct a phylogeny that is best explained by your input data.

We can now start to reconstruct a maximum-likelihood tree from the alignments you received. Make sure you are in the same directory as your trimal file.

`iqtree -s $name.fasta`

`-s` is an option to specify the name of the alignment. This command is always required in order for IQtree to work. 

You will receive three output files after running IQtree. 

    $name.fasta.iqtree
    
    $name.fasta.treefile
    
    $name.fasta.log
    
.iqtree is the main file that you should look at to see the calculations. It also contains view of the final tree.

.treefile is the evolutionary tree that can be viewed in programs like FigTree and iTOL.

.log is a file of the entire command run.
#### Determing what model is best
IQtree supports a wide range of models for DNA, protein, codons, binary and nonbinary alignments. If you are not sure which model to use for your file you can run ModelFinder to figure out which is best.

`-m MFP` is the option to specify the model name to use during the analysis. The special MFP key word stands for ModelFinder Plus, which tells IQ-TREE to perform ModelFinder and the remaining analysis using the selected model. 

#### Assessing branch support with bootstrap estimate
IQtree presents an ultrafast bootstrap estimate that is faster than the standard procedure and provides relatively unbiased bootstrap support estimates. 

To run UF boot

`iqtree -s $name.fasta -m MFP -B 1000`

`-B` specifies the number of bootstrap replicates where 1000 is the minimum number recommended.

Use the space below to show how you will run IQtree

```
```

