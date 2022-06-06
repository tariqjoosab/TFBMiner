"""
Conducts the program execution.
"""

import time
import sys
import glob
import os

import pandas as pd

import generate_chains
import create_args
import single_gene_prediction
import chain_processor


def main():
    t1 = time.time()
    args = create_args.argument_parser()
    inducer = args.compound
    max_chain_length = args.length
    single_gene_operons = args.single_gene_operons 

    if max_chain_length < 2:
        sys.exit("Error: chains cannot be less than 2 enzymes in length.")
    if max_chain_length > 5:
        i = input("The maximum chain length is recommended to be <= 5 for a manageable runtime. Would you still like to continue? (y/n): ")
        if i == 'y':
            pass
        else:
            sys.exit("Program closed.")
    if not os.path.isdir("genome_files"):
        sys.exit("Error: 'genome_files' folder was not found within the program directory.")
    if not os.listdir("genome_files"):
        sys.exit("Error: 'genome_files' folder is empty.")
        
    genome_files = glob.glob("genome_files\*")
    try:
        genome_assemblies = pd.read_csv("genome_assemblies.csv")
        genome_assemblies.drop(columns=genome_assemblies.columns[0], 
            axis=1, 
            inplace=True)
    except FileNotFoundError:
        sys.exit("Error: 'genome_assemblies.csv' was not found within the program directory.")
    
    if single_gene_operons != "y":
        print("{}Identifying enzymatic chains for '{}' with maximum chain length set to {}...{}".format("\n", inducer, max_chain_length, "\n"))
        chains = generate_chains.optimize_chain_identifications(inducer, max_chain_length)
        total_chains = len(chains)
        print("{}{} unique chains were identified.".format("\n", total_chains))
        chain_processor.process_all_chains(chains, total_chains, inducer, genome_assemblies, genome_files, t1)
    else:
        single_gene_prediction.predict_single_gene_operon_biosensors(inducer, genome_assemblies, genome_files, t1)

