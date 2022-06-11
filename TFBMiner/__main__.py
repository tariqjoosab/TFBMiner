"""
Conducts the software execution.
"""

import time
import sys
import glob
import os

import pandas as pd

from TFBMiner import generate_chains, create_args, single_gene_prediction, chain_processor


_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_path(path):
    return os.path.join(_ROOT, 'data', path)


def main():
    t1 = time.time()
    args = create_args.argument_parser()
    inducer = args.compound
    max_chain_length = args.max_chain_length
    single_gene_operons = args.single_gene_operons 
    genome_files_path = args.genome_files_path
    output_path = args.output_path

    if max_chain_length < 2:
        sys.exit("Error: chains cannot be less than 2 enzymes in length.")
    if max_chain_length > 5:
        i = input("The maximum chain length is recommended to be <= 5 for a manageable runtime. Would you still like to continue? (y/n): ")
        if i == 'y':
            pass
        else:
            sys.exit("Program closed.")

    home_path = os.path.expanduser("~")
    if genome_files_path == "unspecified":
        genome_files_path = os.path.join(home_path, "genome_files")
    if output_path == "unspecified":
        output_path = os.path.join(home_path, "TFB_predictions")

    if not os.path.isdir(genome_files_path):
        sys.exit(f"Error: {genome_files_path} was not found.")
    if not os.listdir(genome_files_path):
        sys.exit(f"Error: {genome_files_path} is empty.")
    else:
        genome_files = glob.glob(genome_files_path+"\*")
        genome_assemblies_path = get_path("genome_assemblies.csv")
    try:
        genome_assemblies = pd.read_csv(genome_assemblies_path)
        genome_assemblies.drop(columns=genome_assemblies.columns[0], 
            axis=1, 
            inplace=True)
    except FileNotFoundError:
        sys.exit(f"Error: {genome_assemblies_path} was not found.")
    
    if single_gene_operons != "y":
        print("{}Identifying enzymatic chains for {} with maximum chain length set to {}...{}".format("\n", inducer, max_chain_length, "\n"))
        chains = generate_chains.optimize_chain_identifications(inducer, max_chain_length)
        total_chains = len(chains)
        print("{}{} unique chains were identified.".format("\n", total_chains))
        chain_processor.process_all_chains(chains, total_chains, inducer, genome_assemblies, genome_files, t1, output_path)
    else:
        single_gene_prediction.predict_single_gene_operon_biosensors(inducer, genome_assemblies, genome_files, t1, output_path)


if __name__ == "__main__":
    main()
