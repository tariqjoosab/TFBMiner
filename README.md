# TFBMiner
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

A data acquisition and analysis pipeline for the rapid identification of putative transcription factor-based biosensors.
## Synopsis
TFBMiner predicts putative transcription factor-based biosensors (TFBs) for a compound of interest firstly by identifying enzymes that sequentially catabolize the compound and linking them to form chains. Each chain is then processed to identify whether the enzymes are encoded by any catabolic operons within bacterial genomes, and putative transcriptional regulators of the catabolic operons are predicted and scored based upon a conceptual model of how TFBs are frequently genetically organized. TFBMiner also has an option for predicting TFBs that regulate single genes, rather than genes encoding enzymatic chains.

## Usage
```sh
py -m TFBMiner [-h] [-l L] [-s S] [-g G] [-o O] compound
```

## Options
`compound`: The KEGG COMPOUND database ID of the compound to predict TFBs for.

`-l, --length`: Specify the maximum integer length of the enzymatic chains to generate. Default = 3

`-s, --single_gene_operons`: Specify whether to predict biosensors for genes that encode single enzyme metabolizers, rather than genes encoding enzymatic chains (y/n). Default = n

`-g, --genome_files_path`: Specify the absolute path of the `genome_files` directory, in which the downloaded feature table bacterial genomes will be located. Otherwise, TFBMiner will try to locate it within the user's home directory.

`-o, --output_path`: Specify the absolute path of the desired output directory. Otherwise, TFBMiner will output the results to the user's home directory.

`-h, --help`: Display the software usage, description, options, and guidance in the terminal.

## Examples
```sh 
py -m TFBMiner C00259 -l 3
```
* Predicts TFBs for l-arabinose (ID: `C00259`)
* Identified chains will be up to `3` enzymes in length
* Nothing specified for `-g`, so TFBMiner will default to using the user's home path to find the `genome_files` directory
* Nothing specified for `-o`, so TFBMiner will default to using the user's home path and as a place to store the predictions

```sh
py -m TFBMiner C01494 -l 5 -g C:\Users\user\Desktop\genome_files -o C:\Users\user\Documents\Results
```
* Predicts TFBs for ferulic acid (ID: `C01494`) 
* Identified chains will be up to `5` enzymes in length
* Feature table genomes will be accessed via `C:\Users\user\Desktop\genome_files` (Windows OS)
* Predictions will be output to `C:\Users\user\Documents\Results` (Windows OS)

```sh
py -m TFBMiner C00180 -s y -o /Users/user/Desktop/Results
```
* Predicts TFBs for benzoate (ID: `C00180`)
* Genes that encode single enzyme metabolizers of benzoate will be used to predict TFBs
* Nothing specified for `-g`, so TFBMiner will default to using the user's home path to find the `genome_files` directory.
* Predictions will be output to `/Users/user/Desktop/Results` (Mac OS X)

## Dependencies
* [numpy](https://numpy.org/) (version 1.21.5)
* [pandas](https://pandas.pydata.org/) (version 1.3.5)
* [tqdm](https://github.com/tqdm/tqdm) (version 4.62.3)

## Setup
To install and run TFBMiner, [Python](https://www.python.org/) (version compatibility: >=3.8, <3.11) must first be installed on the user's system. TFBMiner can be installed using the Python package manager ([pip](https://pypi.org/project/pip/)) via the terminal:
```sh
py -m pip install git+https://github.com/tariqjoosab/tfb-miner.git
```

To process identified enzymatic chains, TFBMiner needs access to complete and fully annotated [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)  feature table genomes for all bacteria held on the [KEGG GENOME](https://www.genome.jp/kegg/genome/) database. These genomes can be downloaded from [this Dropbox folder](https://www.dropbox.com/sh/ezo6ahj033cev8b/AADm-bC728rD0l9PTgPA9bgpa?dl=0). Once downloaded, the (unzipped) folder can be placed within the user's home directory (`C:\Users\user` on Windows OS, for instance), which is where TFBMiner will default to searching within to find the folder; this is conducted in an OS-independent manner. Alternatively, one can place the folder within a different directory and specify its absolute path to TFBMiner via the `-g` command-line argument.

While not necessary, one may wish to use up-to-date versions of the bacterial feature table genomes. If so, they can be downloaded from the [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) database in bulk. However, without advanced specification, this will also result in the acquisition of genomes that are not held on KEGG, and using more genomes than necessary may cause a slight performance deficit. To obtain only relevant genomes, one can paste the contents of `search_phrase.txt` into the advanced search builder of NCBI Assembly. These contents consist of each GenBank assembly code of bacteria on KEGG GENOME separated by the `OR` search operator, which therefore specifies to NCBI Assembly to only retrieve these genomes.

## How it works
### Enzymatic chain identification
TFBMiner initially receives the [KEGG COMPOUND](https://www.genome.jp/kegg/compound/) database ID of a compound of interest (`C1`) and uses the KEGG REST API to retrieve the reactions that it's involved in. Subsequently, reactions that metabolize `C1` are identified, and the IDs of their products are used to identify reactions that metabolize their products. Enzymes that catalyse the initial reactions are linked to enzymes that catalyse the subsequent reactions to form chains that sequentially metabolise `C1`. These initial chains are extended by continuing this process until the chains reach the maximum chain length, which is set by the user. Chains of lower lengths that precede the extended chains are not discarded; all identified chains are sent off to the next stage for processing.

### TFB prediction
Each enzymatic chain is processed to identify putative transcriptional regulators of `C1` degradation. This begins by determining whether any genes that encode enzymes within the chain have genetic organisations that are characteristic of catabolic operons. The KEGG REST API is used to retrieve genes that encode each enzyme within the chain and the organisms that possess them, and the results are filtered to leave only organisms that possess all of the enzymes. For each organism, the software uses an internal database to identify the GenBank accession code of its genome, and then searches locally for a feature table genome that contains this accession code in its filename. The genome is then parsed, and the contents are used to predict operons that facilitate `C1` degradation by evaluating the genetic organisations of the relevant genes. If the genes are clustered on the same DNA strand, they are marked as an operon. Putative transcriptional regulators of the operons are predicted by identifying the nearest upstream transcription factor gene on the opposite DNA strand, as this is a highly frequent genetic organisation of TFBs. Each prediction is scored; regulators that are directly upstream of their operons receive the highest score (0), and points are deducted based upon linear distance from the operon and the strand orientations of the genes that are situated in-between. Predictions are ranked in order of their scores and output to `.csv` files.

## Author
Tariq Joosab.

## Acknowledgements
Research supervisors: Dr Erik Hanko & Prof Rainer Breitling.
