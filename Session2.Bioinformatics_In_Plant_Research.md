

###### tags: `target_capture_workshop`
This section of the target capture workshop defines and elaborates on 'bioinformatics' as a term. We will also go through the roles of a bioinformaticist, and what projects a bioinformaticist might work on. The workshop leaders will share what bioinformatics projects they have worked on. Then, we will learn about sequence files and further, learn how to quality check them before doing analysis.
# SESSION 2: BIOINFORMATICS IN PLANT RESEARCH AND SEQUENCE FILES

[TOC]
## 2.1 Steps of Bioinformatic Research

While the bioinformatics portion of a research project is challenging by itself, running programs and writing code typically follows many months of work of other processes. The road to publication in BioInformatics/Genomics is oftentimes long, with many stages, all equally important to the quality of the results.

1. Field Work/Collection
2. Sample Processing
3. Wet Lab
4. Sequencing
5. Dry Lab/ Data Analysis
6. Writing and Publication

Fieldwork and collections takes a long time to plan and perform, especially if your collections will require extensive travel, funding, and permissions. Careful consideration into species selection, permit applications, and storage will be required. As a bioinformaticist, why would this be critical to your role on the project? You should make yourself aware of the condition of your collection- are you using degraded herbarium specimens? Are your samples ancient? Does your species have a reference genome? Do you have enough of a sample to be able to rerun if coverage is low?

Botanists are encouraged to create vouchers for collections made, in order to store and preserve the physical and chemical properties of the samples long-term. Each sample taken during field work must be recorded and given a unique identifier, in order to connect the specimen to the sequencing data. Your work in bioinformatics may have to be retracable to specific specimens.

Wet lab procedures, including enzymatic fragmentation, pooling, concentration, will affect either the quality of data, or the way your data can be interpreted bioinformatically.

Sequencing will be completed by either a third-party, or within your institution. It will be the job of the bioinformaticist to recieve and interpret files when they arrive from the sequencing service.

Dry Lab/ Data Analysis is where a bioinformaticist processes, inteprets, and displays data. This workshop teaches skills needed to complete this step of the research process.

Writing and responding to journal editors may require additional bioinformatic inquiry and tweaks to graphs and charts. As a bioinformaticist, you should assist in revision of matierials and writing.

## 2.2 Bioinformatic Research in Biology
> Bioinformatics is the study of biology using the tools of computer science. Bioinformatics enables the study of proteins, genes and genomes using computer algorithms and computer databases. According to a more general National Institutes of Health definition (see PDF below), bioinformatics is "research, development, or application of computational tools and approaches for expanding the use of biological, medical, behavioral or health data, including those to acquire, store, organize, analyze, or visualize such data." The related discipline of computational biology is "the development and application of data-analytical and theoretical methods, mathematical modeling and computational simulation techniques to the study of biological, behavioral, and social systems."

SOURCE: [Kennedy Krieger Institute](https://www.kennedykrieger.org/research/centers-labs-cores/bioinformatics)

Bioinformatics is inherintly associated to biology research. However, the use of bioinformatics is not limited to the subjects covered in this tutorial. While it has prevalent use in genomics and genetics projects, computational skills utilized in this workshop can be applied to research in nearly all discliplines within biology. Not only is processing data from wet lab part of bioinformatics, but the visualization of data is a major aspect of bioinformatics work. 

## 2.3 Examples of Bioinformatic Projects
Bioinformatics is a merger of Biology, mathematics, and computer sciences. Bioinformatics projects could include nearly any project revolving around the whole genome, scaled all the way down to specific genes and RNA. Bioinformatics can even be used to analyze the structure of a molecule. Bioinformatic projects also include methods development, like creating primers and optimizing data processing.

Instructors at the workshop will tell you about their projects using bioinformatics, and what impacts their bioinformatics work has had on the field of Biology.

## 2.4 Short Reads Sequencing

When we sequence short reads, we mean that the genomic code (DNA) has been fragmented into small chunks, called 'reads'. We do this because we are unable to sequence very long stretches of DNA, and could potentially lose data. In the wetlab procedure, DNA extraction materials are fragmented with an enzyme called fragmentase. The seqments are amplified using PCR, and then paired together in the library preparation step. This is why data rendered from sequencing will have two files for one sample. One is the forward read, and one is the reverse read. 

![SOURCE: https://thesequencingcenter.com/knowledge-base/what-are-paired-end-reads/](https://i.imgur.com/fR9x2vr.jpg)
SOURCE: [The Sequencing Center](https://thesequencingcenter.com/knowledge-base/what-are-paired-end-reads/)

## 2.5 Opening and Viewing a FASTQ and FASTA files

Last week, we downloaded a file from the KEW Tree of Life, which should be in your home directory. In this session we will open it in the command line.

A FASTA file is a file format that stores sequence data on nucleic acids or protein sequences. This file format allows the user to add comments or notations. The FASTA file format is in text format, which means we are able to open and view it.

A FASTQ file is also a text format, but contains the raw or modified output of DNA sequencers, such as the Illumina system. Each sequence in a FASTQ file is known as a "read" and is represented by four lines:

1. The identifier for the read
2. The DNA sequence
3. A placeholder line, usually just a `+` character
4. [Quality scores represented by a single character.](https://en.wikipedia.org/wiki/Phred_quality_score)

Return to the [Kew Tree of Life website](https://treeoflife.kew.org/) and select a species to download the raw sequences. For this exercise it will be best to select a seqeunce where the first column says `PAFTOL`.

After clicking on your chosen species, follow the `Data Access link` to the data repository where the raw sequence data is stored. The files will be named with the extension `.fastq.gz` You will download **two files** for your chosen species - these will be explained in more detail in the next section. 

Right-click on the link to the `.fastq.gz` file and use `wget` to download the file to the `TC_Workshop` directory  you made in the last session. Repeat this for the other file.

**Expected Results**

After downloading both read files for your sample, your `TCWorkshop` directory should look something like this:

```
[cpu-23-37 TC_Workshop]$ ls
INSDC.ERR5034798.Afrofittonia_silvestris.a353.fasta  Test1
SRR7451100_1.fastq.gz                                Test2
SRR7451100_2.fastq.gz
```

'gz' indicates that the file is 'gzipped' or compressed. To view the KEW Tree of life sequence , we will need to unzip the FASTQ file first using the `gunzip` command, followed by the file name. 

Unzipped FASTQ files look like this: 
`ERR5033410_1.fastq` `ERR5033410_1.fastq`

Zipped FASTA files look like this:
`KEWsequence.fasta.gz`

Unzipped FASTA files look like this:
`KEWsequence.fasta`



### Action 2.5.1: gunzip
*Use the `gunzip` command to unzip your files*
```

```

After the file is unzipped, we can view it using the `less` command, followed by the file name. Use the enter key to travel down the file line by line. 

You can also use the commands `head` and `tail`. These commands are convenient for large files, as they print the first 10 lines on the screen, from either the top or the bottom, depending on which command you use.

### Action 2.5.2: Exploring FASTA/FASTQ files
*Use the `head`, `tail`, and `less` command to view the file*
```

```

What do you see in the file? What is the difference between a FASTQ and a FASTA file?

## 2.6 Managing Programs with Conda

It can be a challenge to maintain programs on a remote system, especially if the programs are only distributed as source code and need to be compiled on each system separately. In Bioinformatics, often different workflows require different versions of the same software, and these can also be very difficult to maintain.

You have already seen one solution to this issue (Singularity), and this week we will introduce another solution: `conda`. [From the Conda website](https://docs.conda.io/en/latest/):

> Conda is an open source package management system and environment management system that runs on Windows, macOS and Linux. Conda quickly installs, runs and updates packages and their dependencies. Conda easily creates, saves, loads and switches between environments on your local computer. It was created for Python programs, but it can package and distribute software for any language.

### Installing Conda
To use Conda, you will need to run the basic installation script. We will be using a version of conda known as "Miniconda" that has a small number of packages included. In a web browser navigate to the [MiniConda installation page](https://docs.conda.io/en/latest/miniconda.html). 
 

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

When the installation is finished, you will need to logout from JupyterLab and log back in (because the script modified your `.bashrc` file).

To finish installing `mamba`:

```
conda install mamba
```

### Installing required packages

Conda can install packages from a variety of sources - there are default packages and additional packages in *repositories*. When you install packages with `conda` you need to tell `conda` where to look. For bioinformatics packages, Bioconda is a common respository. You need to set up conda to look in each of these repositories. **Run these commands in order**

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```


### Creating a Conda environment

Mamba is already installed on Cyverse, so you are ready to create your first environment. Each environment will contain a separate version of Python and any other software you choose to install. 

This session, we will be working with the quality filtering programs, so we will name our environment accordingly:

```
conda create -n filter
```


activate it with:

```
conda activate filter
```

Note that your command prompt has changed with `(filter)` at the beginning to indicate what environment you are in. If you ever need to get out of the envirnoment, use `mamba deactivate`.

### Installing Software in a Conda Environment

Make sure you have the HybPiper environment active by checking that your command prompt begeins with `(filter)` Then install `fastp` and `fastqc` with this command:

```
conda install -c bioconda fastp fastqc
```

Here, the `-c` is telling `mamba` to look for a package called `fastp` and `fastqc`. 


## 2.6 Using FastQC to Evaluate Quality of Reads
FastQC is a tool used to perform quality control functions on untrimmed sequence data. It will also provide a set of summary statistics that allow us to identify problems in the data before continuing with the pipeline. FastQC creates a file which can be opened outside of your terminal. You can read more about the function of FastQC on the [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) site. Here are links to [bad](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) and [good](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) outputs.

FastQC accepts multiple input fastp files at once. This is a good opportunity to utilize the 'wildcard' `*` syntax. Using a star in place of a filename, or part of a filename, will tell the ternimal that you want to run all samples that match the remaining descriptors. For example, `*.fasta` will refer to all files that have the extention .fasta. Additionally, SampleNumber*.fasta will refer to all files that start with 'SampleNumber' and end with the .fasta extension. Using a wildcard will be very helpful when you have a large dataset.

Below are the arguments:

>  fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
           [-c contaminant file] seqfile1 .. seqfileN


### Action 2.6.1
*Try writing a line of code below that prompts FastQC to create an output file.*
```


```

On the GUI to the left of the screen, navigate into your current directory and left click on the fastqc output files and download them. You may need to unzip them to view.
## 2.7 Flags and Helper Documents

If you find yourself stuck or confused while writing a line of code while using a program, try using the helper printout. For example, FASTQC has a helper document that explains arguments needed for running the program. Access it using `fastqc -h`.

Some programs may use `man program` or `program --help` for their help documents. It is also worth looking for the program's webpage to view recent documentation -- the developer may also have set up a help forum or e-mail list to get assistance.

Other good sources of information for UNIX systems is [StackOverflow](https://stackoverflow.com/). More specific help for bioinformatics can be found using [biostars](https://www.biostars.org/), [seqanswers](http://seqanswers.com/), or the [bioinformatics reddit forum](https://www.reddit.com/r/bioinformatics/).

When in doubt, it is frequently useful to use a Google search with the name of the program and the error message. Often, these searches will take you directly to a post on one of the sites mentioned above!


## 2.8 Trimming Data
We will use FastP to trim down reads before pairing them. We trim reads for many reasons, such as removing the adapter sequences from a read (adapter sequences are critical to the library prep stage, but cause problems in bioinformatics), or to trim out the low quality, short, sections of reads. We want to keep only high quality and longer reads to input into the pipeline. 

> fastp supports both single-end (SE) and paired-end (PE) input/output.
>* for SE data, you only have to specify read1 input by -i or --in1, and specify read1 output by -o or --out1.
>* for PE data, you should also specify read2 input by -I or --in2, and specify read2 output by -O or --out2.
>* if you don't specify the output file names, no output files will be written, but the QC will still be done for both data before and after filtering.
>* the output will be gzip-compressed if its file name ends with .gz

`fastp -i [untrimmed.forward].fq -I [untrimmed.reverse].fq -o [trimmedoutput.forward].fq -O [trimmedoutput.reverse].fq`

SOURCE: [OpenGene/ FastP Github Page](https://github.com/OpenGene/fastp)

In addition to the simple usage of FastP, we also need to add additional arguments for the type of data we have. The instructors will inform you of additional arguments needed based on the data set you use.


### Action 2.9.1
*Try adapting the simple usage code for your dataset.*
```



```



## 2.9 Comparing trimmed and untrimmed output

Previously, we ran FastQC on the KEW Tree of life data BEFORE we trimmed the reads.

### Action 2.9.1
*Try using FastQC on the trimmed reads*
```

```

Do you see differences between the two outputs?