ConkatSeq Manual -- still in progress
================

#### ConkatSeq is a python program for concatenating amplicon sequencing data back to its chromosomal location in the genome. (Niv describe) 

Table of Contents
-----------------

- [Data Input/Output](#inputandoutput)
- [Installation](#installation)
  1. [Dependencies](##hardware)
  2. [Install](#install)
        *  [ Git](#git)
        *  [Conda](#conda)
        *  [ConkatSeq conda environment](#conkatseq)
        *  [Download ConkatSeq](#downloadconkatseq)     
- [Usage](#usage)
  1. [Generate clustering table](#table) (build_clustering_table.py)
  2. [Filter clustering table](#polish) (filter_clustering_table.py)
  3. [Compute and Graph clustering table](#graph) (conkat_seq.py)

- [Example](#example)
  

## <a name="inputandoutput"></a> Data Input/Output


ConkatSeq uses three python scripts **(build_clustering_table.py, filter_clustering_table.py, conkat_seq.py)** to process **demultiplexed amplicon sequencing data (FASTA format)** and produce a **network graph of highly associated coccurence sequence clusters (GraphML format)**.  


|**Script**|**Input**|**Description**|**Output**|**Description**|
|---|---|---|---|---|
|**build_clustering_table.py**|**1) FASTA sequences (Illumina amplicon sequencing data file(s))**| 1) The amplicon sequencing data file(s) must be demultiplexed or split into its corresponding barcoded numbers of well plate. <br/><br/> *Ex: If the sequencing was performed on a 384 well plate, the demultiplexed data will consist of 384 individual FASTA files representing each barcoded sample well in the plate.*|**1) (sample_name)_OTU.txt**  <br/><br/> **2) (sample_name)_OTU.fna** | 1)Domain clustering table <br/><br/>  2) Domain centroid sequences
|**filter_clustering_table.py**| 
|**conkat_seq.py**|



## <a name="installation"></a> Installation

### <a name="hardware"></a> i. Dependencies

ConkatSeq currently supports any machine with Linux operating system that has Python 2.7 and will require users to use a conda environment for installing the necessary software dependencies for ConkatSeq to run. 

|**Hardware**|**Software**|**Package Manager**|
|---|---|---|
|**1) [Unix](https://www.linux.org/pages/download/)** <br/><br/> **(Any Mac or Linux OS machine)**|**1) [python 2.7](https://www.python.org/download/releases/2.7/)**  <br/><br/>  **2) [biopython](https://biopython.org/)** <br/><br/>  **3) [pandas](https://pandas.pydata.org)** <br/><br/> **4) [scipy](https://www.scipy.org/)** <br/><br/>  **5) [matplotlib](https://matplotlib.org/)** <br/><br/> **6) [statsmodel](https://www.statsmodels.org/stable/index.html)** <br/><br/> **7) [networkx](https://networkx.github.io/)** <br/><br/> **8) [vsearch](https://github.com/torognes/vsearch)**  | **1) [conda](https://conda.io/en/latest/)** 

### <a name="install"></a> ii. Install

All installations will take place on the command line via the terminal.

#### <a name="git"></a>       **a) Git:** 

All related material required to run ConkatSeq can be found in this github repository which can be downloaded via git. 

The following link provides instructions to install git on a Mac or Linux OS: **https://git-scm.com/downloads**

#### <a name="conda"></a>     **b) Conda:** 

Download and install the conda installer for Python 2.7 from the Anaconda distribution manager:

**[Mac](https://www.anaconda.com/distribution/#macos)**: 

```
terminal$ wget  https://repo.anaconda.com/archive/Anaconda2-2019.03-MacOSX-x86_64.sh
terminal$ bash  Anaconda2-2019.03-MacOSX-x86_64.sh
```

**[Linux](https://www.anaconda.com/distribution/#linux)**: 

```
terminal$ wget  https://repo.anaconda.com/archive/Anaconda2-2019.03-Linux-x86_64.sh
terminal$ bash  Anaconda2-2019.03-Linux-x86_64.sh
```
**NOTE:** Make sure you are in the path where the Anaconda2-2019.03-******-x86_64.sh file is located and with actual name of the file.

#### <a name="conkatseq"></a>  <a name="conkatseq2"></a>    **c) ConkatSeq conda environment:** 

Create a conda environment with Python 2.7 dedicated to run Conkatseq:

```
terminal$ conda create --name conkatseq python=2.7
```
Actiavte the ConkatSeq conda environment:

```
terminal$ conda activate conkatseq
(conkatseq)[terminal]$
```

**NOTE:** To deactivate or get out of the ConkatSeq conda environment, type the following: ```(conkatseq)[terminal]$ conda deactivate conkatseq```

Install the required software dependencies inside the ConkatSeq conda environment: 

```
(conkatseq)[terminal]$ conda install -c anaconda pandas networkx statsmodels scipy
(conkatseq)[terminal]$ conda install -c conda-forge biopython matplotlib  
(conkatseq)[terminal]$ conda install -c bioconda vsearch
```

#### <a name="downloadconkatseq"></a> <a name="download2"></a> **d) Download ConkatSeq:** 

Use the git command to obtain this ConkatSeq repository which will contain all the necessary material required for ConkarSeq to run:

```
(conkatseq)[terminal]$ git clone https://github.com/brady-lab-rockefeller/conkat_seq
(conkatseq)[terminal]$ ls conkat_seq
conkat_seq  README.md  resources
(conkatseq)[terminal]$ cd conkat_seq
(conkatseq)[terminal conkat_seq]$ ls
build_clustering_table.py  conkat_utils.py   filter_clustering_table.py  helpers.pyc  ref
conkat_seq.py              conkat_utils.pyc  helpers.py                  __init__.py
```
## <a name="usage"></a> Usage

ConkatSeq consist of three main components divided in three scripts which is how the tool runs. The following will describe how each script is being used with its data and options.

**NOTE:** To run ConkatSeq you must activate its conda environment created in the [ConkatSeq conda environment](#conkatseq2) section earlier and be inside the conkat_seq folder which was [downloaded from this repository](#download2): 

### <a name="table"></a>  <a name="build2"></a>  **i) Generate clustering table** 

#### a) build_clustering_table.py:

The build_clustering_table.py is the first step and component in the ConkatSeq tool. 

The following command shows how to run the build_clustering.py script with its required options:

```
(conkatseq)[terminal conkat_seq]$ python build_clustering_table.py  -i INPATH -o OUTPATH -s SAMPLE_NAME -l STRIP_LEFT -t TRUNCATE -c CLUSTER_ID
```

|**Required Script options**|**Option use**|
|---|---|
|**-i INPATH**|**ENTIRE** folder path of the input file(s) must be provided|
|**-o OUTPUT**|**ENTIRE** folder path of where the output files should be located must be provided|  
|**-s SAMPLE_NAME**|name to provide for the outputs|  
|**-l STRIP_LEFT**|numerical value that represent the number if bases to remove at the left side of the read (usually primer length) |  
|**-t TRUNCATE**|numerical value that represent the number of bases to keep for each read|  
|**-c CLUSTER_ID**|numerical value in decimal point (ex: 0.95) that represents sequence holomolgy identity threshold for clustering|  


For more information about the options type the following that will provide descriptions for every existing parameter.:
```
(conkatseq)[terminal conkat_seq]$ python build_clustering_table.py -h
usage: build_clustering_table.py [-h] -i INPATH -o OUTPATH -s SAMPLE_NAME -l
                                 STRIP_LEFT -t TRUNCATE -c CLUSTER_ID
                                 [--host_path HOST_PATH] [--threads THREADS]
                                 [--verbose] [--remove_files]

build_clustering_table script

optional arguments:
  -h, --help            show this help message and exit
  -i INPATH, --inpath INPATH
                        subpool demultiplexed read files directory
  -o OUTPATH, --outpath OUTPATH
                        output directory for processed files
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        sample name for output files
  -l STRIP_LEFT, --strip_left STRIP_LEFT
                        number of bases to strip (typically primer length)
  -t TRUNCATE, --truncate TRUNCATE
                        number of bases to keep
  -c CLUSTER_ID, --cluster_id CLUSTER_ID
                        identity threshold for amplicon clustering, default =
                        0.95
  --host_path HOST_PATH
                        path to host refrence fasta file if read filtering is
                        required, default = False
  --threads THREADS     number of threads to use by vsearch, default = 1
  --verbose             increase output verbosity
  --remove_files        do not keep processed read files




### <a name="polish"></a> **ii) Filter clustering table** 

### <a name="graph"></a> **iii) Compute and Graph clustering table** 
```
#### b) filter_clustering_table.py:

The following command shows how to run the filter_clustering_table.py script which uses the outputs produced from the first component from the [build_clustering_table.py](#build2)  with its required options:

```
python filter_clustering_table.py  -i INPATH -s SAMPLE_NAME -mrs MIN_READ_SIZE -rst RELATIVE_SIZE_THRESHOLD -msp MIN_SUBPOOLS 
```
|**Required Script options**|**Option use**|
|---|---|
|**-i INPATH**|**ENTIRE** folder path of the input file(s) must be provided|
|**-s SAMPLE_NAME**|name to provide for the outputs|  
|**-mrs MIN_READ_SIZE**|numerical value that represents the number of amplicon reads (ex: 3) to comsider |  
|**-rst  RELATIVE_SIZE_THRESHOLD**|numerical value that represents the pvalue threshold (ex: 0.05) used to remove insignificant reads|  
|**-msp MIN_SUBPOOLS**|numerical value in decimal point (ex: 3) that represents sequence holomolgy identity threshold for clustering|  

```

(conkatseq)[terminal conkat_seq]$ python filter_clustering_table.py -h
usage: filter_clustering_table.py [-h] -i INPATH -s SAMPLE_NAME -mrs
                                  MIN_READ_SIZE -rst RELATIVE_SIZE_THRESHOLD
                                  -msp MIN_SUBPOOLS [--threads THREADS]
                                  [--verbose]

build_clustering_table script

optional arguments:
  -h, --help            show this help message and exit
  -i INPATH, --inpath INPATH
                        directory containing the domain clustering table
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        sample name for output files
  -mrs MIN_READ_SIZE, --min_read_size MIN_READ_SIZE
                        only cosndier amplicon with more reads than
                        min_read_size
  -rst RELATIVE_SIZE_THRESHOLD, --relative_size_threshold RELATIVE_SIZE_THRESHOLD
                        relative size threshold for removing amplicons with
                        low reads within clusters, default = 0.05
  -msp MIN_SUBPOOLS, --min_subpools MIN_SUBPOOLS
                        only consider amplicons detected in more than
                        min_subpools subpools
  --threads THREADS     number of threads to use by vsearch, default = 1
  --verbose             increase output verbosity
  ```

#### c)conkat_seq.py:



## <a name="example"></a> Example


https://rockefeller.app.box.com/s/rhrgw13ux6qdns0vax5i4cgyvgtfdyd7




