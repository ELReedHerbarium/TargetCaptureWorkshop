###### tags: `target_capture_workshop`

# SESSION 5: Species Trees

---

[TOC]

---

NOTE: This session expects you to have finished session 4 on gene trees. However, if you do not have gene trees or alignments for several genes, they will be included in the Github repository in the `data` folder:

* `workshop.trimal.tar.gz` : aligned seqeunces for 353 genes and 50 samples
* `workshop_iqtree.tar.gz` : 353 gene trees with bootstrap support

## 5.1 Why do we need a species tree?
Reconstructing species tree is an essential step into assessing relationships among our target species. When we generated trees from our targeted genes, we may notice that each gene tree only represents its own evolutionary history which differs from expected species tree (because of gene duplication, lost, convention, historical hybridization, incomplete lineage sorting, and presence of deep coalescence). In other word, these gene trees may differ from each other and not agree with species tree which reflects the real phylogeny of associated species. In order to acquire a consensus phylogeny to study in related species, we would like to infer species tree using specific methods.

![Figure1DeepCoalescence](https://github.com/gudusanjiao/targetseqphylo/blob/main/miscellaneous/pic2.png?raw=true "Fig 1 Deep Coalescence")

Figure 1. Deep coalescence may cause discordance between gene tree and species tree.

More information about why gene trees often do not agree with species tree, you can refer to [Szöllősi et al.](https://doi.org/10.1093/sysbio/syu048)
## 5.2 How to infer phylogenies of species?
There are two major ways in inferring species tree. First and simple one is to use a concatenated dataset, which refers to **concatenation**. To apply this method, every gene will be added together and aligned into a supermatrix. Then, further phylogenetic inferences will be deployed onto this concatenated gene which will be treated as a single gene complex. 

However, evidence implied that the concatenation is not always correctly reconstruct the phylogeny. It lacks modeling deep coalescence which creates disagreement with the true phylogeny. This is mainly caused by a majority of genes may not share the same history as species do. Thus, when using concatenation method, signals from discordant genes would outweighing the less coalescence ones.

In order to solve the problem of deep coalescence in species tree reconstruction, several methods have been developed to either **model the process of deep coalescence** or taking **shortcuts** where deep coalescence trees still consist with real species tree. Full modeling of deep coalescence is more accurate but cannot take large datasets due to the complexity of computational time. Since genomic level sequencing is more and more prevalent nowadays (including our target sequencing strategy), shortcuts methods are becoming popular because they are usually effective and relatively reliable.

![image](https://user-images.githubusercontent.com/16470742/143940257-bec4517c-f3df-4b03-8e38-fc6239745896.png "Fig 2 Species Tree Methods")

Figure 2. Schematic of the concatenation and coalescent paradigms in phylogenetics. ([Liu et al.](https://nyaspubs.onlinelibrary.wiley.com/doi/full/10.1111/nyas.12747))

## 5.3 Practicing species tree method using ASTRAL
In today's workshop, we will use [ASTRAL](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) first, one of the shortcuts methods, to infer our species tree from our gene tree datasets. 

ASTRAL is a Java phylogenetics program that uses the frequency of quartets in gene trees to estimate a species tree. It offers a one command line solution under Linux. Comparing these quartets is not computational intense so this program may finish in a short time. However, to be consistent to other sessions in our workshop, we will still run this under HPCC job submitting system.

Load up the Cyverse Discovery environment and access the terminal. As before we will create a `mamba` environment for this session and install packages.

```
mamba create -n species_trees
mamba activate species_trees
```

In some cases you may need to switch to environments used in previous sessions; for this you can use `mamba deactivate.`

We will install packages using `mamba` as in previous sessions:


```bash
mamba install astral-tree newick_utils msaconverter
```

This provides the executable script `astral` which is a wrapper around the Java program; [Newick Utils](https://github.com/tjunier/newick_utils), a set of helper scripts for working with phylogenetic tree files; and [msaconverter](https://github.com/linzhi2013/msaconverter), a handy script for converting DNA sequence alignment files.


### Collapsing Gene Trees

The ASTRAL input file is a file containing gene trees inferred from the same taxa. This may be a file you completed in a previous section, or you can download a set of gene tree files from the GitHub repository for this workshop.

First, create a single file containing all the gene tree files:


```
cat *.treefile > genetrees.tre
```


The gene trees you download (or produce from IQTree in the previous session) are fully bifurcating. However, ASTRAL will consider every branch when estimating its Maximum Quartet Species Tree. There is some debate about whether unsupported branches in gene trees can contribute to poorly resolved ASTRAL trees.

We can use the `nw_ed` tool in Newick Utils to collapse gene tree branches below a certain gene tree bootstrap value threshold. For example, to collapse all branches below 50% bootstrap support:

```
nw_ed input.tre 'i & b<50' o > output.tre
```

### Action 5.3.1

Use `nw_ed` to collapse all of the gene trees branches that have less than 33% support and save this as a new file.

```

```

Next, we will run ASTRAL. View the help options with `astral -h`. At a minimum we need to provide the input file of gene trees and the location of the output file species tree:


```bash
astral -i genetrees.tre -o astral.tre
```

This program should be finished within one minute. If ASTRAL is taking longer than one minute, inspect the log on the screen. One common error is that the names of samples are not consistent across genes. Output file will also be in a Newick tree-form text. You may download it to local for further viewing and editing.



### Action 5.3.2
*Use ASTRAL to generate results from gene trees*
```

```

## 5.4 Practicing species tree method using SVDQuartets (on PAUP)
SVDQuartets is also a method assessing quartets trees in input file. Instead of using gene trees, it takes single nucleotides onto quartets to calculate SVD (Singular Value Decomposition). Thus, the input file is a concatenated FASTA (actually should be coverted into NEXUS format, an universal phylogeny prep format for a lot of programs). To clarify, although the input is concatenated but this is not a concatenation method.

![image](https://github.com/gudusanjiao/targetseqphylo/blob/main/miscellaneous/svdq.png?raw=true)

Figure 3. Quartets and SVDQuartets. ([Source](https://www.asc.ohio-state.edu/kubatko.2//SVDQ_Intro.pdf))

PAUP (Phylogenetic Analysis Using Parsimony [and other methods]) is a software platform integrating different types of phylogenetic analysis, including SVDQuartets. Although PAUP has a graphical interface (GUI) version, it is currently only available for Windows, so we will need to download a command line version for Linux.

### Action 5.4.1: Downloading PAUP
Download [PAUP](http://phylosolutions.com/paup-test/) for CentOS Linux from the website using `wget`. 

Once you have downloaded PAUP you will need to unzip the file with `gunzip` and then make PAUP executable using `chmod`:

`chmod +x paup4a168_centos64`


### Action 5.4.2: Concatenating alignments

PAUP requires the sequences in a different file format called NEXUS. First, we can take all of the alignments from an earlier session and merge them together, inserting gaps where there are missing sequences at certain genes. Navigate to the directory that has your alignment files from the earlier session. 

There is a script included in HybPiper to help with making a concateated "supermatrix" but it can be hard to access if HybPiper was installed with `conda`. The script can be downloaded directly:

```
wget https://raw.githubusercontent.com/mossmatters/HybPiper/master/hybpiper/fasta_merge.py
```

Run the `fasta_merge.py` script with Python. Assuming you are in a directory containing a lot of fasta files to merge:

```
python fasta_merge.py --fastafiles *.fasta > supermatrix.fasta
```


Convert your merged FASTA file to NEXUS using `msaconverter`:

```
msaconverter -i supermatrix.fasta -o supermatrix.nexus -p fasta -q nexus
```


### Action 5.4.3: Running PAUP

PAUP uses its own command line - access it by executing `paup`. If you are in the same directory as the `paup` executable, you will need to specify the path - for example `./paup`

Your command prompt will change to be `paup>`

Load your NEXUS file into memory using the `execute` command in `paup`:

```
execute supermatrix.nexus
```

We use the `svdq` command in `paup` to calculate the SVDQuartets tree. First, display the options for `svdq`

```
help svdq
```

The parameters we want to modify are:
- `bootstrap` (to calculate a bootstrap support using 100 pseudoreplicates)
- `nthreads` (to make use of the multiple CPUs on your Cyverse instance)
- `treemodel` (to set whether we assume sites evolve under the multispecies coalescent)

```
svdq boostrap=standard nthreads=8 treemodel=shared
```

This will generate a tree similar to using a **concatenation method** by assuming all sites evolved under the same tree. 

After program finishing running, the bootstrap consensus tree will be displayed. This is the summary of support across the 100 pseudoreplicates. We can save this tree to a file using the `conTree` command:

```
conTree strict=no majRule=yes percent=10 treeFile=svdq.shared.tre
```

### Action 5.4.4
In this section, repeat the SVDQ analysis, but use `treemodel=mscoalescent`  This will let program to assess coalescense using SVDQuartets. Save the output tree as before.

## 5.5 Interpreting Species Trees results
After generating species tree from ASTRAL and PAUP, we may view and edit it through some tree viewer software by rerooting, collapsing clades, re-coloring, enhancing branches, showing branch length etc. The goal of customizing it is to better convey our result to readers. [FigTree](https://github.com/rambaut/figtree/releases) is a Java-based tree viewer offering numbers of options in customize phylogenetic trees. If you want to install it on your own computer, you need a [Java](https://www.java.com/en/download/) environment installed first.

When you load Newick-form Tree file (.tre) in FigTree, you will be asked for the label of node/branch values. Notice that, support values in ASTRAL species tree means the **Local Posterior Probability**, which indicates how well the shown quartet is compared to the other two possible quartets at that branch. The branch lengths are **coalescent units** in related to ancestral population size and the amount of change found on the gene trees. 

### Action 5.5.1
*Download Figtree. Unzip it. And transfer ASTRAL output to local folder to view it with FigTree*

As we said above, FigTree offers tons of options to modify our initial input file. There are several buttons on the top of the window for quick functions to modify a tree. First action to try is to reroot our output species tree based on known outgroup. To apply this, you may select the outgroup branch to make it highlighted, then click "Reroot" button on the top. You may notice that the order of our species tree seems to be changed. Actually, the topology of this tree is not altered at all, which means the relationships between each species on the tree remain the same.

### Action 5.5.2
*Reroot species tree by selected outgroup*

After rerooting our tree, we also want to show the node statistics for each bifurcations. In FigTree, we can click on triangle beside "Node Labels" on sidebar to expand some options available for node information display. In order to display probability on node, we will choose our named label (software asked when you opened the tree file) on "Display". Then, we can customize the font, size, color and format etc. to better convey the information.

### Action 5.5.3
*Display node values in a nice-looking way*

For some phylogeny trees in papers, they often highlighted some branches to inform readers for specific reasons such as shared ancestry, co-evolution events, and homoplasy traits. These highlightings usually enhance the effectiveness in conveying information. To practise this, we can click on triangle beside "Tip Labels" on side bar. Then, we select one branch and click on color option on the side. Choose a color that you think will be appropriate for your tree.

### Action 5.5.4
*Highlight branches that belong to the same taxanomy groups*

There are a lot of other options for you to customize your phylogenetic trees in FigTree. We would like to leave this oppotunity of exploring it to you. You can try but not be limited to change the size of tip labels, change size of branches, change appearance of the tree into cladogram, make highlighted blocks, add scales to tree, and rename tips etc. The ultimate goal is to help yourself interpret the tree and help other better understand what you want to exhibit.

### Action 5.5.5
*Customize both your ASTRAL and SVDQuatets trees in FigTree and try to interpret and compare*

