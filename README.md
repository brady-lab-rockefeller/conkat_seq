Welcome to the ConkatSeq Documentation Manual
=============================================

#### ConkatSeq is a python program for concatenating amplicon sequencing data back to its chromosomal location in the genome. 

Table of Contents
-----------------

- [Description](#about)
  - [Concept](#concept)
  - [Algorithm](#algorithm)
- [Input](#input)
- [Output](#output)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [3 modules](##)
  - [](##)
- [Repeat graph](#graph)
- [benchmark](#performance)
- [Algorithm Description](#algorithm)

## <a name="about"></a> Description

### <a name="concept"></a> Concept

### <a name="algorithm"></a> Algorithm

## <a name="input"></a> Input

Illumina amplicon sequencing data file(s) in FASTQ format.

### NOTE: The amplicon sequencing data file(s) must be demultiplexed or split into its corresponding barcoded number of well plate. 

#### Ex: If the sequencing was performed on a 384 well plate, the demultiplexed data will consist of 384 individual FASTQ files representing each barcoded sample well in the plate.

## <a name="input"></a> Output

Network graph of highly associated coccurence clusters in GraphML format. 



