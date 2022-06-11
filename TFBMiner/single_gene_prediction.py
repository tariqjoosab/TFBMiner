"""
Uses data acquisition and biosensor prediction algorithms
to predict biosensors that regulate single genes.
"""

from itertools import chain
import multiprocessing
import concurrent.futures
import time
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

from TFBMiner import acquire_data, biosensor_predictor, output


def predict_single_gene_operon_biosensors(compound, genome_assemblies, genome_files, t1, output_path):
    """
    Identifies enzymes that metabolize a compound to find possible 
    single-gene operons, and predicts and outputs potential biosensors.
    """
    enzymes = []
    reactions = acquire_data.identify_reactions(compound)
    num_reactions = len(reactions)
    cores = multiprocessing.cpu_count()
    print("{}Identifying enzymes that metabolize '{}'...{}".format('\n', compound, '\n'))
    
    processes = cores if num_reactions>=cores else 2 if num_reactions>=2 else 1 if num_reactions==1 else None
    if processes is not None:
        if processes >= 2:
            # Identifies catabolic enzymes concurrently.
            reactions = np.array(reactions, dtype=object)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                futures = [executor.submit(acquire_data.identify_single_metabolizers, data, compound) for data in np.array_split(reactions, processes)]            
                enzymes = [future.result() for future in futures]
        elif processes == 1:
            enzymes = acquire_data.identify_single_metabolizers(reactions[0])
    else:
        enzymes = []

    # Removes repeated enzymes.
    enzymes = list(chain(*enzymes))
    enzymes = list(dict.fromkeys(enzymes))
    num_enzymes = len(enzymes)
    print("{}{} unique enzymes were identified as metabolizers of {}.".format('\n', num_enzymes, compound))
    if num_enzymes == 0:
        sys.exit()

    print("{}Processing {} enzymes...{}".format('\n', num_enzymes, '\n'))
    total_biosensors = 0
    for n in tqdm(range(num_enzymes)):
        enzyme = enzymes[n].lower()
        encoders = acquire_data.retrieve_encoders(enzyme)
        if encoders is not None:
            encoders_df = pd.DataFrame(encoders, 
            columns=["Organism", str(enzyme) + "_gene(s)"])
            if (not encoders_df.empty) and (len(encoders_df.columns) > 1):
                df_cols = encoders_df.columns.values.tolist()
                biosensors = biosensor_predictor.optimize_biosensor_predictions(encoders_df, genome_assemblies, genome_files, single_gene_operons=True)
                num_biosensors = len(biosensors)
                if num_biosensors > 0:
                    total_biosensors += num_biosensors
                    output.output_predictions(biosensors, compound, df_cols, output_path)
    t2 = time.time()
    if total_biosensors > 0:
        print("{}Processing is complete. {} potential biosensors were identified for {}. Results have been deposited to {}. Total runtime: {}s.".format("\n", total_biosensors, compound, output_path, round(t2-t1, 2)))
    else:
        print("{}Processing is complete. {} potential biosensors were identified for {}. Total runtime: {}".format("\n", total_biosensors, compound, round(t2-t1, 2)))
