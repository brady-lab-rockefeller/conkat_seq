ConkatSeq Manual
================

#### ConkatSeq is a python program for concatenating amplicon sequencing data back to its chromosomal location in the genome. (Niv describe) 

Table of Contents
-----------------

- [Data Input/Output](#inputandoutput)
- [Installation](#installation)
  1. [Dependencies](##hardware)
  2. [Install](#install)
        *  [ Git:](#git)
        *  [Conda:](#conda)
        *  [ConkatSeq:](#conkatseq)
        *  [Install ConkatSeq:](#installconkatseq)     
- [Usage](#usage)
  1. [Generate clustering table](#table)
  2. [Polish clustering table](#polish)
  3. [Compute and Graph clustering table](#graph)

- [Example](#example)
  

## <a name="inputandoutput"></a> Data Input/Output


ConkatSeq uses three python scripts (build_clustering_table.py, filter_clustering_table.py, conkat_seq.py) to process demultiplexed amplicon sequencing data (FASTA format) and produce a network graph of highly associated coccurence sequence clusters (GraphML format).  


|**Script**|**Input**|**Description**|**Output**|**Description**|
|---|---|---|---|---|
|**build_clustering_table.py**|**1) FASTA sequences (Illumina amplicon sequencing data file(s))**| 1) The amplicon sequencing data file(s) must be demultiplexed or split into its corresponding barcoded numbers of well plate. <br/><br/> *Ex: If the sequencing was performed on a 384 well plate, the demultiplexed data will consist of 384 individual FASTA files representing each barcoded sample well in the plate.*|**1) (sample_name)_OTU.txt**  <br/><br/> **2) (sample_name)_OTU.fna** | **1)Domain clustering table** <br/><br/>  **2) Domain centroid sequences** 
|**filter_clustering_table.py**| 
|**conkat_seq.py**|



## <a name="input"></a> Installation

### <a name="hardware"></a> i. Dependencies

ConkatSeq currently supports any machine with Linux operating system and will require users to use a conda environment for installing the necessary software dependencies for ConkatSeq to run. 

|**Hardware**|**Software**|**Package Manager**|
|---|---|---|
|**1) [Unix](https://www.linux.org/pages/download/)** <br/><br/> **(Any Mac or Linux OS machine)**|**1) [python 2.7](https://www.python.org/download/releases/2.7/)**  <br/><br/>  **2) [biopython](https://biopython.org/)** <br/><br/>  **3) [pandas](https://pandas.pydata.org)** <br/><br/> **4) [scipy](https://www.scipy.org/)** <br/><br/>  **5) [matplotlib](https://matplotlib.org/)** <br/><br/> **6) [statsmodel](https://www.statsmodels.org/stable/index.html)** <br/><br/> **7) [networkx](https://networkx.github.io/)** <br/><br/> **8) [vsearch](https://github.com/torognes/vsearch)**  | **1) [conda](https://conda.io/en/latest/)** 

### <a name="install"></a> ii. Install

All installations will take place on the command line via the terminal.

#### <a name="git"></a>**a) Git:** 

All related scripts required to run ConkatSeq can be found in this github repository which can be obtained via git. 

The following links provides instructions to install git on a Mac or Linux OS: https://git-scm.com/downloads

#### <a name="conda"></a>**b) Conda:** 

**Download** and **install** the conda installer for Python 2.7 from the Anaconda distribution manager:

Mac (https://www.anaconda.com/distribution/#macos): 

```
terminal$ wget  https://repo.anaconda.com/archive/Anaconda2-2019.03-MacOSX-x86_64.sh
terminal$ bash  Anaconda2-2019.03-MacOSX-x86_64.sh
```

Linux (https://www.anaconda.com/distribution/#linux): 

```
terminal$ wget  https://repo.anaconda.com/archive/Anaconda2-2019.03-Linux-x86_64.sh
terminal$ bash  Anaconda2-2019.03-Linux-x86_64.sh
```
**NOTE:** Make sure you are in the path where the Anaconda2-2019.03-******-x86_64.sh file is located and with actual name of the file.

**3) ConkatSeq conda environment:** 

Create a conda environment with Python 2.7 dedicated for Conkatseq to run:

```
terminal$ conda create --name conkatseq python=2.7
```
Actiavte the ConkatSeq conda environment:

```
terminal$ conda activate conkatseq
terminal(conkatseq)$
```

Install the required software dependencies inside the ConkatSeq conda environment: 

```
terminal(conkatseq)$ conda install -c anaconda pandas networkx statsmodels scipy
terminal(conkatseq)$ conda install -c conda-forge biopython matplotlib  
terminal(conkatseq)$ conda install -c bioconda vsearch
```

**4) Install ConkatSeq:** 

Use the git command to obtain this ConkatSeq 

```
terminal(conkatseq)$ conda install -c anaconda pandas networkx statsmodels scipy
terminal(conkatseq)$ conda install -c conda-forge biopython matplotlib  
```

## <a name="input"></a> Usage



## <a name="input"></a> Example


https://rockefeller.app.box.com/s/rhrgw13ux6qdns0vax5i4cgyvgtfdyd7




