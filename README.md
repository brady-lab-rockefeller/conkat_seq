ConkatSeq Manual
================

#### ConkatSeq is a python program for concatenating amplicon sequencing data back to its chromosomal location in the genome. (Niv describe) 

Table of Contents
-----------------

- [Input](#input)
- [Output](#output)
- [Installation](#installation)
  - [Hardware dependencies](#hardware)
  - [Software dependencies](#software)
  - [Install](#install)
  
- [Usage](#usage)
  1. [Generate clustering table](#table)
  2. [Polish clustering table](#polish)
  3. [Compute and Graph clustering table](#graph)

- [Example](#example)
  

## <a name="input"></a> Data Input/Output

|**Script**|**Input**|**Description**|**Output**|**Description**|
|---|---|---|---|---|
|**build_clustering_table.py**|**1) FASTA sequences (Illumina amplicon sequencing data file(s))**| 1) The amplicon sequencing data file(s) must be demultiplexed or split into its corresponding barcoded numbers of well plate. <br/><br/> *Ex: If the sequencing was performed on a 384 well plate, the demultiplexed data will consist of 384 individual FASTA files representing each barcoded sample well in the plate.*|**1) (sample_name)_OTU.txt**  <br/><br/> **2) (sample_name)_OTU.fna** | **1)Domain clustering table ** 2) Domain centroid sequences
|**filter_clustering_table.py**| 
|**conkat_seq.py**|


## <a name="output"></a> Output

Network graph of highly associated coccurence clusters in GraphML format. 

## <a name="input"></a> Installation

### <a name="hardware"></a> Hardware Dependencies

Currently ConkatSeq works in the following operating systems:

  -MacOS
  -Unix
  

### <a name="hardware"></a> Software Dependencies

## <a name="input"></a> Usage

## <a name="input"></a> Example


https://rockefeller.app.box.com/s/rhrgw13ux6qdns0vax5i4cgyvgtfdyd7




