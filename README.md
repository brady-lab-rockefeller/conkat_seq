ConkatSeq Manual
================

#### ConkatSeq is a python program for concatenating amplicon sequencing data back to its chromosomal location in the genome. (Niv describe) 

Table of Contents
-----------------

- [Data Input/Output](#inputandoutput)
- [Installation](#installation)
  - [Dependencies](#hardware)
  - [Install](#install)
  
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

### <a name="hardware"></a> Dependencies


|**Hardware**|**Software**|**Package Manager**|
|---|---|---|
|**Linux** <br/><br/> **(Any Mac or Unix OS machine)**|**[python 2.7]**(https://www.ubuntu.com/)  <br/><br/>  **biopython** <br/><br/>  **pandas** <br/><br/> **scipy** <br/><br/>  **matplotlib** <br/><br/> **statsmodel** <br/><br/> **networkx** <br/><br/> **vsearch**  | **conda** 

### <a name="hardware"></a> Software Dependencies

## <a name="input"></a> Usage

## <a name="input"></a> Example


https://rockefeller.app.box.com/s/rhrgw13ux6qdns0vax5i4cgyvgtfdyd7




