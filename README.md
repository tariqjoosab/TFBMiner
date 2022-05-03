# TFBMiner
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A data acquisition and analysis pipeline for the rapid identification of putative transcription factor-based biosensors.
## Description

### Synopsis
TFBMiner predicts putative transcription factor-based biosensors (TFBs) for a compound of interest firstly by identifying enzymes that sequentially catabolize the compound and linking them to form chains. Each chain is then processed to identify whether the enzymes are encoded by a catabolic operon, and putative transcriptional regulators of the catabolic operons are identified and scored based upon a conceptual model of how TFBs are organized within bacterial genomes.

### Enzymatic chain identification
TFBMiner initially receives the KEGG COMPOUND database ID of a compound of interest (denoted C1) and uses the KEGG REST API to retrieve data regarding which reactions C1 is involved in. Subsequently, reactions that catabolize the compound are identified, and the IDs of each product (C2) are used to identify reactions that catabolize C2. Enzymes that catalyse the initial reactions are linked to enzymes that catalyse the subsequent reactions to form chains that sequentially processes C1. This process continues until chains reach the maximum chain length, which is set by the user. Each chain will be output to the console, keeping the user updated during this stage.

### TFB prediction
Each enzymatic chain is processed to identify putative transcriptional regulators of C1 degradation. This begins by determining whether any genes that encode enzymes within the chain have genetic organisations that are characteristic of a catabolic operon. The KEGG REST API is used to retrieve lists of genes that encode each enzyme within the chain and the organisms that possess them, and this data is filtered to leave only organisms that possess all of the enzymes. 
For each organism, the program uses a local database file (genome_assemblies.csv) in the program directory to identify the GenBank accession code of its genome, and then searches within a subdirectory (genome_files) for a feature table genome that contains this accession code in its filename. The genome is then read, and the data are used to predict operons that facilitate C1 degradation and their putative transcriptional regulators based upon a conceptual model of their relative organisational characteristics within genomes. Each prediction is then scored based upon how well its organisation fits the model. 
Data for each chain is then output to a csv file within a directory named chainlength=x, where x is the length of the chain, and each of these directories will be within a parent directory named C1_results. Predictions will be ranked in order of their scores, with a score of 0 being the highest achievable score. A progress bar in the console will keep the user updated on the progress of this stage.

## Usage
```sh
TFBMiner [-h] [compound] [-l length]
```

## Options
```-h, --help``` 

Display the program usage, description, options, and guidance in the terminal.

```compound```

Enter the KEGG COMPOUND database ID of the compound to perform computations for.

```-l, --length```

Specify the maximum length of the enzymatic chains. This value can currently range between 2 and 4.

## Examples

```sh 
TFBMiner.py C00259 -l 3
```

Performs computations for l-arabinose. Chains of length 2 and 3 will be identified and processed.

```sh
TFBMiner.py C01494 -l 4
```

Performs computations for ferulic acid. Chains of length 2, 3, and 4 will be identified and processed.

## Requirements
- numpy (version: 1.21.5)
- pandas (version: 1.3.5)
- tqdm (version: 4.62.3)

## Setup
To process identified enzymatic chains, complete and fully annotated GenBank feature table genomes of all bacterial genomes held on the KEGG GENOME database should be downloaded and placed within a folder named genome_files. These genomes can be downloaded from the following Dropbox folder: https://www.dropbox.com/sh/ezo6ahj033cev8b/AADm-bC728rD0l9PTgPA9bgpa?dl=0. 

However, it is recommended to use up-to-date versions. Bacterial feature table genomes can be downloaded from the NCBI Assembly database in bulk. However, doing this will also result in the acquisition of genomes that are not held on KEGG, and using more genomes than necessary may cause a slight performance deficit. To obtain only relevant genomes, one can paste the contents of search_phrase.txt into the advanced search builder of NCBI Assembly. These contents consist of each GenBank assembly code of bacteria on KEGG GENOME separated by the “OR” search operator, which therefore specifies to NCBI Assembly to only retrieve these genomes. This search phrase was created on 28/02/2022, and therefore only includes KEGG Genomes that were hosted prior to this date. If one wishes to retrieve genomes added after this date, their GenBank accession codes can simply be added to this search phrase.

## Author
Tariq Joosab.
