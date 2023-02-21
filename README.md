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

## Installation
These are beginner's instructions for using TFBMiner. 

### Dependencies
* Anaconda/miniconda
download minicoda if you dont have it already
https://docs.conda.io/en/latest/miniconda.html

### Windows
Open the Anaconda/miniconda prompt

### Linux
Open a terminal

### Create a virtual environment
Create a virtual environment using the following command  (for this part it doesn't matter where you are in your file system). If it asks if you want to proceed press y.
> conda create --name tfbMiner python=3.10

Activate the environment, you should see "(base)" on the left of the conda prompt be replaced with "(tfbMiner)". It this does not happen then software that you install could affect the result of your computer.
> conda activate tfbMiner

Install the requirements with the following commands 
> conda install pip
> conda install tqdm=4.62.3
> conda install numpy=1.21.5
> conda install pandas=1.5.2

### Download the code/data
Install the code by going to the github page and clicking Code, Download Zip. Unzip the zip-file in the directory that you want to work in https://github.com/RuthStoney/TFBMiner - this will be changed to https://github.com/UoMMIB/TFBMiner ones everything's totally finalized

Download the data by going to the Dropbox page and clicking download. This contains a lot of files so make a folder called genome_files in the directory that you want to work in and extract the data into there.
https://www.dropbox.com/sh/ezo6ahj033cev8b/AADm-bC728rD0l9PTgPA9bgpa?dl=0

### Final steps befor running the code
Make a results folder in the directory that you want to work in.

Within the conda prompt you need to us the "cd" (change directory) command to get to the correct place within your file system.
> cd (dont press enter yet!)

Open a windows explorer and drag the TFBMiner-main file into the prompt, it should look like this
> (tfbMiner) C:\Users>cd C:\Users\ruths\workCode\TFBMiner-main (press enter)

Double check that you are in the TFBMiner-main folder and the tfbMiner environment. The left of your terminal should look something like this:
> (tfbMiner) C:\Users\ruths\workCode\TFBMiner-main>

### Run the code! 
Use the -g flag to indictate the genome data and the -o file to indicate the results folder that you made
> python -m TFBMiner C00259 -g C:\Users\ruths\workCode\TFBMiner-main\genome_files\ -o C:\Users\ruths\workCode\TFBMiner-main\results\

### Trouble shooting
It didn't work - make sure that you are in the correct location within your file system by checking the file path on the left of the Anaconda prompt. The last folder should be TFBMiner

It returned 0 potential biosensors - this suggests that it couldn't find the genome_files. Did you put a slash at the end of your file paths? Double check the filepath provided and that you are running the code from the TFBMiner-main folder (see previous help point).


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


## Setup - details
To process identified enzymatic chains, TFBMiner needs access to complete and fully annotated [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)  feature table genomes for all bacteria held on the [KEGG GENOME](https://www.genome.jp/kegg/genome/) database. These genomes can be downloaded from [this Dropbox folder](https://www.dropbox.com/sh/ezo6ahj033cev8b/AADm-bC728rD0l9PTgPA9bgpa?dl=0). Once downloaded, the (unzipped) folder can be placed within the user's home directory (`C:\Users\user` on Windows OS, for instance), which is where TFBMiner will default to searching within to find the folder; this is conducted in an OS-independent manner. Alternatively, one can place the folder within a different directory and specify its absolute path to TFBMiner via the `-g` command-line argument.

While not necessary, one may wish to use up-to-date versions of the bacterial feature table genomes. If so, they can be downloaded from the [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) database in bulk. However, without advanced specification, this will also result in the acquisition of genomes that are not held on KEGG, and using more genomes than necessary may cause a slight performance deficit. To obtain only relevant genomes, one can paste the contents of `search_phrase.txt` into the advanced search builder of NCBI Assembly. These contents consist of each GenBank assembly code of bacteria on KEGG GENOME separated by the `OR` search operator, which therefore specifies to NCBI Assembly to only retrieve these genomes.

## How it works
### Enzymatic chain identification
TFBMiner initially receives the [KEGG COMPOUND](https://www.genome.jp/kegg/compound/) database ID of a compound of interest (`C1`) and uses the KEGG REST API to retrieve the reactions that it's involved in. Subsequently, reactions that metabolize `C1` are identified, and the IDs of their products are used to identify reactions that metabolize their products. Enzymes that catalyse the initial reactions are linked to enzymes that catalyse the subsequent reactions to form chains that sequentially metabolise `C1`. These initial chains are extended by continuing this process until the chains reach the maximum chain length, which is set by the user. Chains of lower lengths that precede the extended chains are not discarded; all identified chains are sent off to the next stage for processing.

### TFB prediction
Each enzymatic chain is processed to identify putative transcriptional regulators of `C1` degradation. This begins by determining whether any genes that encode enzymes within the chain have genetic organisations that are characteristic of catabolic operons. The KEGG REST API is used to retrieve genes that encode each enzyme within the chain and the organisms that possess them, and the results are filtered to leave only organisms that possess all of the enzymes. For each organism, the software uses an internal database to identify the GenBank accession code of its genome, and then searches locally for a feature table genome that contains this accession code in its filename. The genome is then parsed, and the contents are used to predict operons that facilitate `C1` degradation by evaluating the genetic organisations of the relevant genes. If the genes are clustered on the same DNA strand, they are marked as an operon. Putative transcriptional regulators of the operons are predicted by identifying the nearest upstream transcription factor gene on the opposite DNA strand, as this is a highly frequent genetic organisation of TFBs. Each prediction is scored; regulators that are directly upstream of their operons receive the highest score (0), and points are deducted based upon linear distance from the operon and the strand orientations of the genes that are situated in-between. Predictions are ranked in order of their scores and output to `.csv` files.

## Author
Tariq Joosab & Ruth Stoney

## Acknowledgements
Research supervisors: Dr Erik Hanko & Prof Rainer Breitling.
