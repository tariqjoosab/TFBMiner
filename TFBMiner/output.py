"""
Formats and outputs the biosensor prediction data.
"""

import os
import re
import csv


def output_predictions(biosensors, inducer, df_cols):
    """
    Formats the predicted biosensors data to be output
    as .csv files in specific directories.
    """
    biosensors.sort(key=lambda x: x.regulator_score, reverse=True)
    num_cols = len(df_cols)

    # Prepares output directory names for the data.
    dir_1 = "Results"
    dir_2 = str(inducer) + "_results"
    if num_cols > 2:
        dir_3 = "chainlength=" + str(num_cols-1)
    else:
        dir_3 = "single-enzyme_predictions"

    # Prepares name for .csv file that will hold predictions
    # for a specific enzymatic chain.
    enzyme_colnames = df_cols[1:]
    last_ec_num_idxs = [colname.rfind(re.match('.+([0-9])[^0-9]*$', colname).group(1)) for colname in enzyme_colnames]    
    filename = inducer + "(" + "-".join("ec" + colname[3:last_ec_num_idxs[enzyme_colnames.index(colname)]+1] for colname in enzyme_colnames) + ").csv"
    
    # Formats the data for .csv files.
    header = ["Organism_code"] + [df_cols[x] for x in range(1, num_cols)] + ["Operon", 
                "Regulator", 
                "Regulator_score", 
                "Regulator_annotation"]
    
    data = [[biosensor.organism_code] + 
            [biosensor.genes[x] for x in range(1, num_cols)] +
            [" ".join(str(gene) for gene in biosensor.operon),
            biosensor.regulator,
            biosensor.regulator_score,
            biosensor.regulator_annotation] 
            for biosensor in biosensors]
    
    # Creates directories and outputs
    # data as .csv files within them.
    try:
        path = os.path.join(dir_1, dir_2, dir_3, filename)
        with open(path, "w", encoding="UTF8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(data)
            
    except FileNotFoundError:
        dirs = os.path.join(dir_1, dir_2, dir_3)
        os.makedirs(dirs)
        path = os.path.join(dir_1, dir_2, dir_3, filename)
        with open(path, "w", encoding="UTF8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(data)