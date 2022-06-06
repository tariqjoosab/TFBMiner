# TFBMiner
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

A data acquisition and analysis pipeline for the rapid identification of putative transcription factor-based biosensors.
## Description

### Synopsis
TFBMiner predicts putative transcription factor-based biosensors (TFBs) for a compound of interest firstly by identifying enzymes that sequentially catabolize the compound and linking them to form chains. Each chain is then processed to identify whether the enzymes are encoded by any catabolic operons, and putative transcriptional regulators of the catabolic operons are predicted and scored based upon a conceptual model of how TFBs are organized within bacterial genomes. TFBMiner also has an option for predicting TFBs that regulate single genes, rather than genes encoding enzymatic chains.

### Enzymatic chain identification
TFBMiner initially receives the [KEGG COMPOUND](https://www.genome.jp/kegg/compound/) database ID of a compound of interest (denoted C1) and uses the KEGG REST API to retrieve data regarding which reactions C1 is involved in. Subsequently, reactions that catabolize the compound are identified, and the IDs of each product (C2) are used to identify reactions that catabolize C2. Enzymes that catalyse the initial reactions are linked to enzymes that catalyse the subsequent reactions to form chains that sequentially processes C1. This process continues until chains reach the maximum chain length, which is set by the user. Each chain will be output to the console, keeping the user updated during this stage.

### TFB prediction
Each enzymatic chain is processed to identify putative transcriptional regulators of C1 degradation. This begins by determining whether any genes that encode enzymes within the chain have genetic organisations that are characteristic of a catabolic operon. The KEGG REST API is used to retrieve genes that encode each enzyme within the chain and the organisms that possess them, and the results are filtered to leave only organisms that possess all of the enzymes. For each organism, the program uses a local database file (genome_assemblies.csv) in the program directory to identify the [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) accession code of its genome, and then searches within a subdirectory (genome_files) for a feature table genome that contains this accession code in its filename. The genome is then read, and the contents are used to predict operons that facilitate C1 degradation and their putative transcriptional regulators based upon a conceptual model of their relative organisational characteristics within genomes. Each prediction is then scored based upon how well its organisation fits the model. The highest achievable score is 0, and points are deducted due to divergence from the model. Predictions for each chain are then output to a csv file within a directory named 'chainlength=x', where x is the length of the chain, and each of these directories will be within a parent directory named 'C1_results'. Predictions will be ranked in order of their scores. A progress bar in the console will keep the user updated on the progress of this stage.

## Usage
```sh
TFBMiner.py [-h] [-l LENGTH] [-s SINGLE_GENE_OPERON] compound
```

## Options
```-h, --help```: Display the program usage, description, options, and guidance in the terminal.

```compound```: Enter the KEGG COMPOUND database ID of the compound to predict TFBs for.

```-l, --length```: Specify the maximum length of the enzymatic chains.

```-s, --single_gene_operons```: Choose whether to predict biosensors for rare, potential single-gene operons (y/n).

## Examples
To predict TFBs for l-arabinose using generated chains up to 3 enzymes in length:
```sh 
TFBMiner.py C00259 -l 3
```

To predict TFBs for ferulic acid using generated chains up to 5 enzymes in length:

```sh
TFBMiner.py C01494 -l 5
```

To predict TFBs for ferulic acid that may regulate single genes:

```sh
TFBMiner.py C01494 -s y
```

## Dependencies
- [numpy](https://numpy.org/) (version 1.21.5)
- [pandas](https://pandas.pydata.org/) (version 1.3.5)
- [tqdm](https://github.com/tqdm/tqdm) (version 4.62.3)

## Setup
To run TFBMiner, its dependencies and [Python 3.10.1](https://www.python.org/downloads/release/python-3101/) must be installed. A requirements.txt file has been provided to allow one to easily install the dependencies using the Python package manager ([pip](https://pypi.org/project/pip/)):

```sh
pip install requirements.txt
```

To process identified enzymatic chains, complete and fully annotated GenBank feature table genomes of all bacterial genomes held on the KEGG GENOME database should be downloaded and placed within a folder named genome_files. These genomes can be downloaded from [this Dropbox folder](https://www.dropbox.com/sh/ezo6ahj033cev8b/AADm-bC728rD0l9PTgPA9bgpa?dl=0). 

However, it is recommended to use up-to-date versions. Bacterial feature table genomes can be downloaded from the [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) database in bulk. Doing this will also result in the acquisition of genomes that are not held on KEGG, and using more genomes than necessary may cause a slight performance deficit. To obtain only relevant genomes, one can paste the contents of search_phrase.txt into the advanced search builder of NCBI Assembly. These contents consist of each GenBank assembly code of bacteria on [KEGG GENOME](https://www.genome.jp/kegg/genome/) separated by the “OR” search operator, which therefore specifies to NCBI Assembly to only retrieve these genomes. This search phrase was created on 28/02/2022, and therefore only includes KEGG Genomes that were hosted prior to this date. If one wishes to retrieve genomes added after this date, their GenBank accession codes can simply be added to this search phrase.

## Author
Tariq Joosab.

## Acknowledgements
Research supervisors: Dr Erik Hanko & Prof Rainer Breitling.
