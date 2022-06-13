"""
Contains functions for applying the biosensor prediction algorithms on 
identified enzymatic chains.
"""

import time

import pandas as pd
from tqdm import tqdm

from TFBMiner import acquire_data, output, biosensor_predictor, identify_metabolizers

class MetabolizerProcessor:
    def __init__(self, inducer, genome_assemblies, genome_files, t1, output_path, metabolizers, total_metabolizers):
        self.inducer = inducer
        self.genome_assemblies = genome_assemblies
        self.genome_files = genome_files
        self.t1 = t1
        self.output_path = output_path
        self.metabolizers = metabolizers
        self.total_metabolizers = total_metabolizers

    def process_chains(self):
        """
        Iterates through a collection of enzymatic chains and processes 
        them to predict potential biosesors for the compound they metabolize.
        """
        print("{}Processing {} chains...{}".format("\n", self.total_metabolizers, "\n"))
        total_biosensors = 0
        for n in tqdm(range(self.total_metabolizers)):
            #num_biosensors = self.process_chain(self.metabolizers[n])
            all_genes = {}
            chain = self.metabolizers[n]
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
                        biosensors = biosensor_predictor.execute_biosensor_predictions(filtered_encoders_df, self.genome_assemblies, self.genome_files)
                        # If biosensors were predicted, they are ranked in order
                        # of their scores and formatted for data output.
                        num_biosensors = len(biosensors)
                        if num_biosensors > 0:
                            output.output_predictions(biosensors, self.inducer, df_cols, self.output_path)
                            total_biosensors+=num_biosensors
                except KeyError:
                    pass
        t2 = time.time()
        if total_biosensors > 0:
            print("{}Processing is complete. {} potential biosensors were identified for {}. Results have been deposited to {}. Total runtime: {}s.".format("\n", total_biosensors, self.inducer, self.output_path, round(t2-self.t1, 2)))
        else:
            print("{}Processing is complete. {} potential biosensors were identified for {}. Total runtime: {}s".format("\n", total_biosensors, self.inducer, round(t2-self.t1, 2)))
            alt = input("Would you like to predict biosensors for potential single-gene operons, instead? Predictions are more likely to be made, but at expense of lower prediction accuracy. (y/n)")
            if alt == "y":
                self.t1 = time.time()
                enzymes, total_enzymes = identify_metabolizers.MetabolizerIdentifier(self.inducer).execute_single_metabolizer_identification()
                self.metabolizers = enzymes
                self.total_metabolizers = total_enzymes
                total_biosensors = self.process_single_metabolizers()

    def process_single_metabolizers(self):
        """
        Iterates through a collection of single enzyme metabolizers and processes 
        them to predict potential biosesors for the compound they metabolize.
        """
        print("{}Processing {} enzymes...{}".format("\n", self.total_metabolizers, "\n"))
        total_biosensors = 0
        for n in tqdm(range(self.total_metabolizers)):
            enzyme = self.metabolizers[n].lower()
            encoders = acquire_data.retrieve_encoders(enzyme)
            if encoders is not None:
                encoders_df = pd.DataFrame(encoders, 
                columns=["Organism", str(enzyme) + "_gene(s)"])
                if (not encoders_df.empty) and (len(encoders_df.columns) > 1):
                    df_cols = encoders_df.columns.values.tolist()
                    biosensors = biosensor_predictor.execute_biosensor_predictions(encoders_df, self.genome_assemblies, self.genome_files, single_gene_operons=True)
                    num_biosensors = len(biosensors)
                    if num_biosensors > 0:
                        total_biosensors += num_biosensors
                        output.output_predictions(biosensors, self.inducer, df_cols, self.output_path)
        t2 = time.time()
        if total_biosensors > 0:
            print("{}Processing is complete. {} potential biosensors were identified for {}. Results have been deposited to {}. Total runtime: {}s.".format("\n", total_biosensors, self.inducer, self.output_path, round(t2-self.t1, 2)))
        else:
            print("{}Processing is complete. {} potential biosensors were identified for {}. Total runtime: {}s".format("\n", total_biosensors, self.inducer, round(t2-self.t1, 2)))