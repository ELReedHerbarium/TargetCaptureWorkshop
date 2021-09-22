# SESSION 2: BIOINFORMATICS IN PLANT RESEARCH AND SEQUENCE FILES

[TOC]
## 2.1 Kew Tree of Life Explorer
Last session, we used the Kew Tree of Life Explorer to download and transfer data to the TTU HPCC. The Kew Tree of Life Explorer is important to genomicists globally, as their data is stored and available for use by other scientists after publication. 

> The Kew Tree of Life Explorer is an output of the Plant and Fungal Trees of Life Project (PAFTOL) at the Royal Botanic Gardens, Kew. PAFTOL aims to discover and disseminate the evolutionary history of all plant and fungal genera. The evolutionary tree of life is fundamental to our understanding of the natural world. Comparative studies of DNA sequence data have revolutionised our knowledge of the tree of life, however many gaps remain. In collaboration with partners from around the world, PAFTOL is addressing this challenge by generating, compiling and analysing genomic data for all ca. 13,900 flowering plant and 8,200 fungal genera to build novel trees of life at unprecedented scale.

SOURCE: [Kew Royal Botanic Gardens](https://treeoflife.kew.org/about)

Currently, the Tree of Life Explorer is focused on mostly flowering plants, with some gymnosperms available. In the future, fungal trees will be added to the repository.
## 2.2 Steps of Bioinformatic Research

While the bioinformatics portion of a research project is challenging by itself, running programs and writing code typically follows many months of work of other processes. The road to publication in BioInformatics/Genomics is oftentimes long, with many stages, all equally important to the quality of the results.

1. Field Work/Collection
2. Sample Processing
3. Wet Lab
4. Sequencing
5. Dry Lab/ Data Analysis
6. Writing and Publication

Fieldwork and collections takes a long time to plan and perform, especially if your collections will require extensive travel and funding. Careful consideration into species selection, permit applications, and storage will be required. As a bioinformaticist, why would this be critical to your role on the project? You should make yourself aware of the condition of your collection- are you using degraded herbarium specimens? Do you have enough of a sample to be able to rerun if recovery is poor?

Botanists are encouraged to create vouchers for collections made, is order to store and preserve the physical and chemical properties of the samples long-term. Each sample taken during field work must be recorded and given a unique identifier, in order to connect the specimen to the sequencing data. Your work in bioinformatics may have to be retracable to specific specimens.

Wet lab procedures, including enzymatic fragmentation, pooling, concentration, will affect either the quality of data, or the way your data can be interpreted bioinformatically.

Sequencing will be completed by either a third-party, or within your lab. It will be the job of the bioinformaticist to recieve and interpret files when they arrive from the sequencing service.

Dry Lab/ Data Analysis is where a bioinformaticist processes, inteprets, and displays data. This workshop teaches skills needed to complete this step of the research process.

Writing and responding to journal editors may require additional bioinformatic inquiry and tweaks to graphs and charts. As a bioinformaticist, you should assist in revision of matierials and writing.

## 2.3 Bioinformatic Research in Biology
> Bioinformatics is the study of biology using the tools of computer science. Bioinformatics enables the study of proteins, genes and genomes using computer algorithms and computer databases. According to a more general National Institutes of Health definition (see PDF below), bioinformatics is "research, development, or application of computational tools and approaches for expanding the use of biological, medical, behavioral or health data, including those to acquire, store, organize, analyze, or visualize such data." The related discipline of computational biology is "the development and application of data-analytical and theoretical methods, mathematical modeling and computational simulation techniques to the study of biological, behavioral, and social systems."

SOURCE: [Kennedy Krieger Institute](https://www.kennedykrieger.org/research/centers-labs-cores/bioinformatics)

Bioinformatics is inherintly associated to biology research. However, the use of bioinformatics is not limited to the subjects covered in this tutorial. While it has prevalent use in genomics and genetics projects, computational skills utilized in this workshop can be applied to research in nearly all discliplines within biology. Not only is processing data from wet lab part of bioinformatics, but the visualization of data is a major aspect of bioinformatics work. 

## 2.4 Examples of Bioinformatic Projects
## 2.5 Opening a FASTA file
A FASTA file is a file format that stores sequence data on nucleic acids or protein sequences. This file format allows the user to add comments or notations. The FASTA file format is in text format, which means we are able to open and view it.

To view the sequences in the test dataset files we have provided, we will need to unzip them first. Unzipped FASTA files look like this:

`
testdataset.fastq
`

While zipped files could look like this

`
testdataset.fq.gz
`

'gz' indicates that the file is 'gzipped'. We can unzip the file using the `gunzip` command, followed by the file name.

```

```

After the file is unzipped, we can view it using the `less` command, followed by the file name. Use the enter key to travel down the file line by line

```

```

What do you see in the file?
## 2.6 Short Reads Sequencing

When we sequence short reads, we mean that the genomic code (DNA) has been fragmented into small chunks, called 'reads'. We do this because we are unable to sequence very long stretches of DNA, and could potentially lose data. In the wetlab procedure, DNA extraction materials are fragmented with an enzyme called fragmentase. The seqments are amplified using PCR, and then paired together in the library preparation step. This is why data rendered from sequencing will have two files for one sample. One is the forward read, and one is the reverse read. 

![SOURCE: https://thesequencingcenter.com/knowledge-base/what-are-paired-end-reads/](https://i.imgur.com/fR9x2vr.jpg)
SOURCE: [The Sequencing Center](https://thesequencingcenter.com/knowledge-base/what-are-paired-end-reads/)

## 2.7 Using FastQC to Evaluate Quality of Reads
FastQC is a tool used to perform quality control functions on untrimmed sequence data. It will also provide a set of summary statistics that allow us to identify problems in the data before continuing with the pipeline. FastQC creates a file which can be opened outside of your terminal. You can read more about the function of FastQC and see example outputs on the [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) site.

FastQC accepts multiple input fastp files at once. This is a good opportunity to utilize the 'wildcard' `*` syntax. Using a star in place of a filename, or part of a filename, will tell the ternimal that you want to run all samples that match the remaining descriptors. For example, `*.fasta` will refer to all files that have the extention .fasta. Additionally, SampleNumber*.fasta will refer to all files that start with 'SampleNumber' and end with the .fasta extension. Using a wildcard will be very helpful when you have a large dataset.

Try writing a line of code below that prompts FastQC to create an output file. Below are the arguments:

>  fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
           [-c contaminant file] seqfile1 .. seqfileN

```


```
## 2.8 Flags and Helper Documents

If you find yourself stuck or confused while writing a line of code while using a program, try using the helper printout. For example, FASTQC has a helper document that explains arguments needed for running the program. Access it using `fastqc -h`.

Helper documents for many functions in command line can be accessed by the -h argument. If in-depth help, find the developer's GitHub page.
## 2.9 Trimming Data
We will use FastP to trim down reads before pairing them. We want to keep only high quality reads to input into the pipeline. 

> fastp supports both single-end (SE) and paired-end (PE) input/output.
>* for SE data, you only have to specify read1 input by -i or --in1, and specify read1 output by -o or --out1.
>* for PE data, you should also specify read2 input by -I or --in2, and specify read2 output by -O or --out2.
>* if you don't specify the output file names, no output files will be written, but the QC will still be done for both data before and after filtering.
>* the output will be gzip-compressed if its file name ends with .gz

`fastp -i [untrimmed.forward].fq -I [untrimmed.reverse].fq -o [trimmedoutput.forward].fq -O [trimmedoutput.reverse].fq`

SOURCE: [OpenGene/ FastP Github Page](https://github.com/OpenGene/fastp)

In addition to the simple usage of FastP, we also need to add additional arguments for the type of data we have. The instructors will inform you of additional arguments needed based on the data set you use.

Try adapting the simple usage code for your dataset.

```



```



