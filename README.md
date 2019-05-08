# CONKAT-seq [work in progress]

CONKAT-seq (co-occurrence network analysis of targeted sequences) is targeted sequencing workflow that enables the exploration of rare biosynthetic gene clusters in complex metagenomes. CONKAT-seq is designed to reconstruct the clustered organization of biosynthetic domains in the metagenome based on the statistical analysis of amplicon co-occurrences in a partitioned library of large insert metagenomic clones. Briefly, high molecular weight DNA from soil samples are extracted and cloned to construct large insert metagenomic library which preserves the linkage between co-clustered genes. Library clones are randomly partitioned into hundreds of wells (subpools), and DNA sequence encoding for biosynthetic domains of interest are amplified using barcoded primers. Amplicon de-barcoding identifies the positioning of each biosynthetic domain within the array of subpools and co-occurrence frequencies of biosynthetic domain variants across all subpools are recorded. Pairwise statistical analysis of domain co-occurrence (Fisher’s exact test) identifies domain pairs that show strong linkage and is used to assign amplicon variants into distinct biosynthetic domain networks. Finally, CONKAT-seq predictions are visualized as networks, where nodes represent sequence variants of the targeted biosynthetic domains and edges link domains predicted to be physically co-clustered in the metagenomic DNA. 

Table of Contents
-----------------
- [How it works](#inputandoutput)
- [Installation & Dependencies](#installation)
- [Usage](#usage)
  1. [Generate clustering table](#table) (build_clustering_table.py)
  2. [Filter clustering table](#polish) (filter_clustering_table.py)
  3. [Compute and Graph clustering table](#graph) (conkat_seq.py)
- [Example](#example)
  
## How it works

CONKAT-seq requires 3 processing steps to process subpool-demultiplexed amplicon sequencing data (FASTA format) to predicted networks of chromosomally co-clustered biosynthetic domains (GraphML format).

#### Pre-processing steps (repeat for every targeted domain amplicon dataset)

**build_clustrering_table.py** 
Pre-processing (primer removal, length trimming, dereplication) of subpool-demultiplexed amplicon reads. Processed reads from all library subpools are clustered using VSEARCH implementation of the USEARCH algorithem. Each cluster contains a set of highly similar amplicon sequences (<95% identity) originating from one or more library subpools.

- input:
    1. Folder containg amplicon sequencing fasta file(s). Files must be demultiplexed according to subpool amplicon barcode. For example, if the targeted domain amplifcation was performed on a 384 subpools library (i.e., 384 PCR reactions), the demultiplexed data will consist of 384 individual FASTA files representing each one subpool sample.
    
- output
    1. Domain amplicons clustering table in a UCLUST-format tabbed text format [sample_name.txt]
    2. Sequences of cluster centroids in a FASTA format	[sample_name.fna]

**filter_clustrering_table.py** 
Parasing of the domain clustering table into a dataframe and filtering of domain varinats with low read counts or low number of subpool occurrences.
- input
    1. Domain amplicons clustering table in a UCLUST-format tabbed text format [sample_name.txt]
    2. Sequences of cluster centroids in a FASTA format	[sample_name.fna]
       
- Output
    1. Filtered domain clustering dataframe [sample_name.csv]

#### Network analysis (Once per metagenomic library. Can integrate multiple domain amplicon datasets.)

**conkat_seq.py**
Pairwise statisical analysis of pairwise domain co-occurances. To identify pairs of biosynthetic domains that originate from physically clustered metagenomic DNA, a 2x2 contingency table (the number of subpools containing both domain variants, one of the two only, or none of them) is constructed for each pair of domain sequence variants the co-occurrence significance is computed using Fisher's exact test. Pairs of domains showing non-random association based on p-value cutoff vlaue are predicted to be physically linked, and hence predicted to belong to the same gene cluster. Based on a pairwise list of statistically significant links a graph representation of domain networks is constructed, where nodes represent cluster of biosynthetic domains and edges link domains that are predicted to be physically co-clustered.

- input
    1. One or more filtered domain clustering dataframe [sample_name.csv]

- output
    1. Predicted networks of chromosomally co-clustered biosynthetic domains in a GraphML format [sample_name.graphml]

## Installation and Dependencies 

CONKAT-seq is available for Linux and MacOS platforms and requires the installation of Python (v2.7.x) and VSEARCH (v2.9.1+). In order to use "clear_host_reads" mode (removal of amplicons matching library host genome, ususally E. coli) BBMap and SAMTOOLS (v3.0.0+) are needed to be in the user path.

#Required Python libraries
- **[biopython](https://biopython.org/)** 
- **[pandas](https://pandas.pydata.org)**
- **[scipy](https://www.scipy.org/)**  
- **[matplotlib](https://matplotlib.org/)** 
- **[statsmodel](https://www.statsmodels.org/stable/index.html)** 
- **[networkx](https://networkx.github.io/)**  

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

#### build_clustering_table.py:

```
usage: build_clustering_table.py [-h] -i INPATH -o OUTPATH -s SAMPLE_NAME -l
                                 STRIP_LEFT -t TRUNCATE -c CLUSTER_ID
                                 [--host_path HOST_PATH] [--threads THREADS]
                                 [--verbose] [--remove_files]
                                 
python build_clustering_table.py  -i INPATH -o OUTPATH -s SAMPLE_NAME -l STRIP_LEFT -t TRUNCATE -c CLUSTER_ID
```
Arguments:  
`-i INPATH` full, absolute path to folder containing the demultiplexed subpool FASTA files

`-o OUTPUT` full, absolute path to where output files will be saved

`-s SAMPLE_NAME` name to provide for the output files

`-l STRIP_LEFT` number of bases to remove at the start of the reads (usually primer length)

`-t TRUNCATE` length of sequence keep following bases removal (shorter sequences are discarded)

`-c CLUSTER_ID`  minimum sequence identity for clustering (a fraction between 0.0 and 1.0, default 0.95)

Optional arguments & flags:  
`--host_path HOST_PATH`  full, absolute path to host reference genome in FASTA format (matching reads are removed)

`--threads THREADS`  threads to be used (default 1)

`--verbose`  increas verbosity

`--remove_files`  remove intermediate processeing files


#### filter_clustering_table.py:

```
usage: filter_clustering_table.py [-h] -i INPATH -s SAMPLE_NAME -mrs
                                  MIN_READ_SIZE -rst RELATIVE_SIZE_THRESHOLD
                                  -msp MIN_SUBPOOLS [--threads THREADS]
                                  [--verbose]

python filter_clustering_table.py  -i INPATH -s SAMPLE_NAME -mrs MIN_READ_SIZE -rst RELATIVE_SIZE_THRESHOLD -msp MIN_SUBPOOLS 
```

Arguments:  
`-i INPATH` full, absolute path to folder containing the input files

`-s SAMPLE_NAME` sample name matching the input files

`-mrs MIN_READ_SIZE` only cosndier amplicon with more reads than min_read_size

`-rst  RELATIVE_SIZE_THRESHOLD` relative size threshold for removing amplicons with low reads within clusters (default 0.05)

`-msp MIN_SUBPOOLS`  only consider amplicons detected in more than min_subpools subpools (default 3)


Optional arguments & flags:  

`--threads THREADS`  threads to be used (default 1)

`--verbose`  increas verbosity


#### conkat_seq.py:

```
usage: conkat_seq.py [-h] -l LIST_OF_CLUSTERING_DATAFRAMES
                     [LIST_OF_CLUSTERING_DATAFRAMES ...] -o OUTPATH
                     [-m MIN_SHARED_OCCURANCES] [-a ALPHA]
                     [--merge_similar_id MERGE_SIMILAR_ID] [--threads THREADS]
                     --flag_edges [--verbose] [--override]

python conkat_seq.py -l LIST_OF_CLUSTERING_DATAFRAMES  -o OUTPATH -a ALPHA -m MIN_SHARED_OCCURANCES  --flag_edges 
```

Arguments:  

`-l LIST_OF_CLUSTERING_DATAFRAMES` list of one or more domain clustering dataframes

`-o OUTPATH` full, absolute path of folder to save results

`-a ALPHA` maximal adjusted p-value threshold (default 10^-6)

`-m MIN_SHARED_OCCURANCES` only analyze domain pairs with co-occurances >min_shared_occurances (default 3)

Optional arguments & flags:

`--merge_similar_id MERGE_SIMILAR_ID` identify threshold for merging similar domains within network (default 0.9)

`--threads THREADS`  threads to be used (default 1)

`--flag_edges` run monte carlo analysis to flag edges potenially affected by index swapping (default False)

`--threads THREADS`  threads to be used (default 1)

`--verbose`  increas verbosity

`--override`  re-write existing files



https://rockefeller.app.box.com/s/rhrgw13ux6qdns0vax5i4cgyvgtfdyd7

## <a name="example"></a> Example

### <a name="initialize"></a> i. Initialize

Make sure you are the conkat_seq folderActivate the conda environment and a

### <a name="runbuild"></a> ii. Run build_clustering_table.py

### <a name="runfilter"></a> iii. Run filter_clustering_table.py

### <a name="runconkat"></a> iv. Run conkat_seq.py







