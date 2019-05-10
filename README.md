<img src="https://github.com/brady-lab-rockefeller/conkat_seq/blob/master/resources/conkat_seq_logo.jpeg" width=60, height=60 align="left" />

# CONKAT-seq 

CONKAT-seq (co-occurrence network analysis of targeted sequences) is targeted sequencing workflow that enables the exploration of rare biosynthetic gene clusters in complex metagenomes. CONKAT-seq is designed to reconstruct the clustered organization of biosynthetic domains in the metagenome based on the statistical analysis of amplicon co-occurrences in a partitioned library of large insert metagenomic clones. Briefly, high molecular weight DNA from soil samples are extracted and cloned to construct large insert metagenomic library which preserves the linkage between co-clustered genes. Library clones are randomly partitioned into hundreds of wells (subpools), and DNA sequence encoding for biosynthetic domains of interest are amplified using barcoded primers. Amplicon de-barcoding identifies the positioning of each biosynthetic domain within the array of subpools and co-occurrence frequencies of biosynthetic domain variants across all subpools are recorded. Pairwise statistical analysis of domain co-occurrence (Fisher’s exact test) identifies domain pairs that show strong linkage and is used to assign amplicon variants into distinct biosynthetic domain networks. Finally, CONKAT-seq predictions are visualized as networks, where nodes represent sequence variants of the targeted biosynthetic domains and edges link domains predicted to be physically co-clustered in the metagenomic DNA. 

Table of Contents
-----------------
- [How it works](#howitworkds)
- [Installation & Dependencies](#Installation)
- [Usage](#usage)
  1. [Generate clustering table](#table) (build_clustering_table.py)
  2. [Filter clustering table](#polish) (filter_clustering_table.py)
  3. [Compute and Graph clustering table](#graph) (conkat_seq.py)
- [Example](#example)
- [Contact](#contact)
  
## <a name="howitworkds"></a> How it works

**CONKAT-seq** requires 3 processing steps to process subpool-demultiplexed amplicon sequencing data (FASTA format) to predicted networks of chromosomally co-clustered biosynthetic domains (GraphML format).

#### Pre-processing steps (repeat for every targeted domain amplicon dataset)
---
**Step 1: build_clustrering_table.py** 

Pre-processing (primer removal, length trimming, dereplication) of subpool-demultiplexed amplicon reads. 

Processed reads from all library subpools are clustered using VSEARCH implementation of the USEARCH algorithm. Each cluster contains a set of highly similar amplicon sequences (<95% identity) originating from one or more library subpools.

- input:
    1. Folder containing amplicon sequencing fasta file(s). Files must be demultiplexed according to subpool amplicon barcode.  
    For example, if the targeted domain amplification was performed on a 384 subpools library (i.e., 384 PCR reactions), the demultiplexed data will consist of 384 individual FASTA files representing each one subpool sample.
    
- output
    1. Domain amplicons clustering table in a UCLUST-format tabbed text format [sample_name.txt]
    2. Sequences of cluster centroids in a FASTA format	[sample_name.fna]

**Step 2: filter_clustrering_table.py** 

Parasing of the domain clustering table into a dataframe and filtering of domain variants with low read counts or low number of subpool occurrences.

- input
    1. Domain amplicons clustering table in a UCLUST-format tabbed text format [sample_name.txt]
    2. Sequences of cluster centroids in a FASTA format	[sample_name.fna]
       
- output
    1. Filtered domain clustering dataframe [sample_name.csv]

#### Network analysis (Once per metagenomic library. Can integrate multiple domain amplicon datasets.)
---
**Steps 3: conkat_seq.py**

Pairwise statistical analysis of pairwise domain co-occurrences. 

To identify pairs of biosynthetic domains that originate from physically clustered metagenomic DNA, a 2x2 contingency table (the number of subpools containing both domain variants, one of the two only, or none of them) is constructed for each pair of domain sequence variants the co-occurrence significance is computed using Fisher's exact test. Pairs of domains showing non-random association based on p-value cutoff value are predicted to be physically linked, and hence predicted to belong to the same gene cluster. Based on a pairwise list of statistically significant links a graph representation of domain networks is constructed, where nodes represent cluster of biosynthetic domains and edges link domains that are predicted to be physically co-clustered.

- input
    1. One or more filtered domain clustering dataframe [sample_name.csv]

- output
    1. Predicted networks of chromosomally co-clustered biosynthetic domains in a GraphML format [sample_name.graphml]

## <a name="Installation"></a> Installation and Dependencies 

CONKAT-seq is available for Linux and MacOS platforms and requires the installation of **[Python (v2.7.x)](https://www.python.org/downloads/release/python-2716/)** and **[VSEARCH (v2.9.1+)](https://github.com/torognes/vsearch)**. 

In order to use "clear_host_reads" mode (removal of amplicons matching library host genome, usually E. coli) **[BBMap](https://sourceforge.net/projects/bbmap/)** and **[SAMTOOLS (v3.0.0+)](http://samtools.sourceforge.net/)** are needed to be in the user path.

## Required Python libraries
- **[biopython](https://biopython.org/)** 
- **[pandas](https://pandas.pydata.org)**
- **[scipy](https://www.scipy.org/)**  
- **[matplotlib](https://matplotlib.org/)** 
- **[statsmodel](https://www.statsmodels.org/stable/index.html)** 
- **[networkx](https://networkx.github.io/)** 

#### **Alternatively users can also install the [conda package manager](https://www.anaconda.com/distribution/) and create a conda environment with all of the required dependencies installed and will be able to use and run CONKAT-seq in the conda environment:**

```
conda install -c anaconda pandas networkx statsmodels scipy
conda install -c conda-forge biopython matplotlib  
conda install -c bioconda vsearch
```

To download CONKAT-seq using Git:
```
git clone https://github.com/brady-lab-rockefeller/conkat_seq
```
## Usage


####  <a name="table"></a> build_clustering_table.py:

```
usage: build_clustering_table.py [-h] -i INPATH -o OUTPATH -s SAMPLE_NAME -l
                                 STRIP_LEFT -t TRUNCATE -c CLUSTER_ID
                                 [--host_path HOST_PATH] [--threads THREADS]
                                 [--verbose] [--remove_files]
                                 
python build_clustering_table.py  -i INPATH -o OUTPATH -s SAMPLE_NAME -l STRIP_LEFT -t TRUNCATE -c CLUSTER_ID
```
**Arguments (required):** 

`-i INPATH` full, **absolute path** to folder containing the demultiplexed subpool FASTA files

`-o OUTPUT` full, **absolute path** to where output files will be saved

`-s SAMPLE_NAME` name to provide for the output files

`-l STRIP_LEFT` number of bases to remove at the start of the reads (usually primer length)

`-t TRUNCATE` length of sequence to keep following bases removal (shorter sequences are discarded)

`-c CLUSTER_ID`  minimum sequence identity for clustering (a fraction between 0.0 and 1.0, default 0.95)

**Optional arguments & flags:**  

`--host_path HOST_PATH`  full, **absolute path** to host reference genome in FASTA format (matching reads are removed)

`--threads THREADS`  threads to be used (default 1)

`--verbose`  increasw verbosity

`--remove_files`  remove intermediate processing files


#### filter_clustering_table.py:

```
usage: filter_clustering_table.py [-h] -i INPATH -s SAMPLE_NAME -mrs
                                  MIN_READ_SIZE -rst RELATIVE_SIZE_THRESHOLD
                                  -msp MIN_SUBPOOLS [--threads THREADS]
                                  [--verbose]

python filter_clustering_table.py  -i INPATH -s SAMPLE_NAME -mrs MIN_READ_SIZE -rst RELATIVE_SIZE_THRESHOLD -msp MIN_SUBPOOLS 
```

**Arguments (required):** 

`-i INPATH` full, **absolute path** to the folder containing the input files

`-s SAMPLE_NAME` sample name matching the input files

`-mrs MIN_READ_SIZE` only consider amplicon with more reads than min_read_size

`-rst  RELATIVE_SIZE_THRESHOLD` relative size threshold for removing amplicons with low reads (fraction of reads in comparison to max within cluster, default 0.05)

`-msp MIN_SUBPOOLS`  only consider amplicons detected in more than min_subpools subpools (default 3)


**Optional arguments & flags:**  
 

`--threads THREADS`  threads to be used (default 1)

`--verbose`  increase verbosity


#### conkat_seq.py:

```
usage: conkat_seq.py [-h] -l LIST_OF_CLUSTERING_DATAFRAMES
                     [LIST_OF_CLUSTERING_DATAFRAMES ...] -o OUTPATH
                     [-m MIN_SHARED_OCCURRENECS] [-a ALPHA]
                     [--merge_similar_id MERGE_SIMILAR_ID] [--threads THREADS]
                     [--flag_edges] [--verbose] [--override]

python conkat_seq.py -l LIST_OF_CLUSTERING_DATAFRAMES  -o OUTPATH -a ALPHA -m MIN_SHARED_OCCURRENECS  --flag_edges 
```

**Arguments (required):** 

`-l LIST_OF_CLUSTERING_DATAFRAMES` list of one or more domain clustering dataframe files

`-o OUTPATH` full, **absolute path**  to save output files 

`-a ALPHA` maximal adjusted p-value threshold (default 10^-6)

`-m MIN_SHARED_OCCURRENECS` only analyze domain pairs with co-occurrences >min_shared_occurrences (default 3)

**Optional arguments & flags:**  

`--merge_similar_id MERGE_SIMILAR_ID` identify threshold for merging similar domains within network (default 0.9)

`--threads THREADS`  threads to be used (default 1)

`--flag_edges` run monte carlo analysis to flag edges potentially affected by index swapping (default False)

`--threads THREADS`  threads to be used (default 1)

`--verbose`  increase verbosity

`--override`  re-write existing files

## <a name="example"></a> Example

This example uses sample data (demultiplexed amplicon sequencing data  (384 fasta files)) made available and can be downloaded from the following link: https://rockefeller.app.box.com/s/rhrgw13ux6qdns0vax5i4cgyvgtfdyd7 being processed in each of 3 processing steps of the CONKAT-seq workflow.

Make sure you are in the conkat_seq folder in order to run:

```
[user@terminal conkat_seq]$ ls
conkat_seq  README.md  resources  
[user@terminal conkat_seq]$ cd conkat_seq/
[user@terminal conkat_seq]$ ls
build_clustering_table.py  conkat_seq.py  conkat_utils.py  conkat_utils.pyc  filter_clustering_table.py  helpers.py  helpers.pyc  __init__.py  ref
```

### <a name="runbuild"></a> i. Run build_clustering_table.py

First, run the build_clustering_table.py script on the sample data. You must provide the **ABSOLUTE** path of the location to the sample data and desired location of the outputs.  

```
[user@terminal conkat_seq]$ python build_clustering_table.py  -i /home/user/conkat_seq/data/ -o /home/user/conkat_seq/output/ -s sample_name -l 23 -t 210 -c 0.95 --threads 20 
Building clustering table from amplicon data...
384 files found...
0 files processed...
100 files processed...
200 files processed...
300 files processed...
Merging de-replicated reads ->  /home/user/conkat_seq/output/sample_name.fna...
Sorting merged reads...
vsearch --threads 64 --sortbylength  /home/user/conkat_seq/output/sample_name.fna --output  /home/user/conkat_seq/output/sample_name_SORTED.fna 
Clustering merged reads...
vsearch --threads 64 --cluster_size  /home/user/conkat_seq/output/sample_name_SORTED.fna --id 0.95 --iddef 1 --sizein --sizeout --centroids  /home/user/conkat_seq/output/sample_name_OTU.fna --uc  /home/user/conkat_seq/output/sample_name_OTU.txt
vsearch v2.13.1_linux_x86_64, 377.3GB RAM, 64 cores
https://github.com/torognes/vsearch

Reading file  /home/user/conkat_seq/output/sample_name_SORTED.fna 100%
146149290 nt in 695949 seqs, min 210, max 210, avg 210
Masking 100%
Sorting by abundance 100%
Counting k-mers 100%
Clustering 100%
Sorting clusters 100%
Writing clusters 100%
Clusters: 69765 Size min 3, max 37910, avg 10.0
Singletons: 34863, 5.0% of seqs, 50.0% of clusters
Amplicon domain clustering table saved ->  /home/user/conkat_seq/output/sample_name_OTU.txt...
Amplicon domain centroid sequences saved ->  /home/user/conkat_seq/output/sample_name_OTU.fna...
```

### <a name="runfilter"></a> ii. Run filter_clustering_table.py

Second, run the filter_clustering_table.py script on the output data produced from the previous step. You must provide the **ABSOLUTE** path of the location to the output data produced from the previous step as the input option for this step. 

```
[user@terminal conkat_seq]$ python filter_clustering_table.py  -i /home/user/conkat_seq/output/  -s sample_name -mrs 3 -rst 0.05 -msp 3 --threads 20
Domain clustering table found -> /home/user/conkat_seq/output/sample_name_OTU.txt...
Domain centroid sequences found -> /home/user/conkat_seq/output/sample_name_OTU.fna
Parsing clustering information...
Populating table features...
Calculating cluster sizes...
69765 clusters found (105315 OTUs)...
Dropping clusters with min_read_size < 3 from table...
25580 clusters found (38537 OTUs)...
Dropping clusters with min_subpools < 3 from table...:
5091 clusters found (16371 OTUs)...
Dropping reads with readSize < max_cluster_readsize*0.05 ...
5091 clusters found (13718 OTUs)...
Remove reads with potenial barcode-swapping between sample plates...
5091 clusters found (13718 rows)...
4940 clusters found (13534 OTUs)...
Parsed and filtered domain clustering dataframe saved -> /home/user/conkat_seq/output/sample_name_OTU.csv...
```

### <a name="runconkat"></a> iii. Run conkat_seq.py

Third, run the conkat_seq.py script on the output data produced from the previous step. You must provide the **ABSOLUTE** path of the location to the domain clustering dataframe(s) (OTU.csv file(s)) produced from the previous step as the input option for this step, and desired location of the outputs.

```
[user@terminal conkat_seq]$ python conkat_seq.py -l /home/user/conkat_seq/output/sample_name_OTU.csv  -o /home/user/conkat_seq/output/ -a 0.05 -m 3  --threads 20
Concatenating domain dataframes...
Calculating domain co-occurrences...
Building domain clustering graph...
Constructing network --> fdr_tsbky 0.05
978 nodes and 1320 edges found...
merging similar domains within networks...
205 networks found... 
0
100
200
Done!
```
### <a name="contact"></a> Contact





