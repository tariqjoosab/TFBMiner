"""
Contains functions for applying the biosensor prediction algorithms on 
identified enzymatic chains.
"""

import time

import pandas as pd
from tqdm import tqdm

from TFBMiner import acquire_data, output, biosensor_predictor, single_gene_prediction


def process_chain(chain, inducer, genome_assemblies, genome_files, output_path):
    """
    Finds organisms that possess all enzymes within a chain and applies 
    the biosensor prediction algorithms to the corresponding genes.
    """
    # Identifies whether there are any genomes that
    # encode all enzymes within the chain.
    all_genes = {}
    for n in range(len(chain)):
        enzyme = chain[n].lower()
        encoders = acquire_data.retrieve_encoders(enzyme)
        if encoders is not None:
            # Organisms and their genes are stored in a dataframe.
            all_genes["d" + str(n)] = pd.DataFrame(encoders, 
            columns=["Organism", str(enzyme) + "_gene(s)"])
    
    num_sets = len(all_genes)
    if num_sets > 1:    
        try:
            # Merges the dataframes to form one that only includes 
            # organisms that possess all enzymes within the chain.
            for x in range(num_sets-1):
                all_genes["d" + str(x+1)] = pd.merge(
                left = all_genes["d" + str(x)], 
                right = all_genes["d" + str(int(x)+1)], 
                left_on = "Organism", 
                right_on = "Organism")

            filtered_encoders_df = all_genes["d" + str(max(range(num_sets)))] 
            if (not filtered_encoders_df.empty) and (len(filtered_encoders_df.columns) > 2):
                df_cols = filtered_encoders_df.columns.values.tolist()
                biosensors = biosensor_predictor.optimize_biosensor_predictions(filtered_encoders_df, genome_assemblies, genome_files)
                # If biosensors were predicted, they are ranked in order
                # of their scores and formatted for data output.
                num_biosensors = len(biosensors)
                if num_biosensors > 0:
                    output.output_predictions(biosensors, inducer, df_cols, output_path)
                    return num_biosensors
        except KeyError:
            pass


def process_all_chains(chains, total_chains, inducer, genome_assemblies, genome_files, t1, output_path):
    """
    Iterates through a collection of enzymatic chains and processes 
    them to predict potential biosesors for the compound they metabolize.
    """
    print("{}Processing {} chains...{}".format("\n", total_chains, "\n"))
    total_biosensors = 0
    for n in tqdm(range(total_chains)):
        num_biosensors = process_chain(chains[n], inducer, genome_assemblies, genome_files, output_path)
        if num_biosensors is not None:
            total_biosensors += num_biosensors
    t2 = time.time()
    if total_biosensors > 0:
        print("{}Processing is complete. {} potential biosensors were identified for {}. Results have been deposited to {}. Total runtime: {}s.".format("\n", total_biosensors, inducer, output_path, round(t2-t1, 2)))
    else:
        print("{}Processing is complete. {} potential biosensors were identified for {}. Total runtime: {}".format("\n", total_biosensors, inducer, round(t2-t1, 2)))
        alt = input("Would you like to predict biosensors for potential single-gene operons, instead? Predictions are more likely to be made, but at expense of lower prediction accuracy. (y/n)")
        if alt == "y":
            t1_new = time.time()
            single_gene_prediction.predict_single_gene_operon_biosensors(inducer, genome_assemblies, genome_files, t1_new, output_path)