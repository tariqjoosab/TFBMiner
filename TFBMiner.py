from itertools import chain
from urllib.error import HTTPError
from urllib.request import urlopen
from tqdm import tqdm
import argparse
import pandas as pd
import numpy as np
import typing as typ
import concurrent.futures
import multiprocessing
import time
import sys
import glob
import os
import csv
import io


def argument_parser():

    """Creates command-line options for the user."""

    parser = argparse.ArgumentParser(description = "TFBMiner: Identifies putative transcription factor-based biosensors for a given compound.")
    parser.add_argument("compound", help="Enter the KEGG Compound ID of the inducer compound.")
    parser.add_argument("-l", "--length", type=int, help="Enter the maximum length of the enzymatic chains.", default=3)

    args = parser.parse_args()

    return args


def get_data(search_term):

    """Retrieves data held in KEGG database entries using KEGG's RESTful API."""

    URL = "http://rest.kegg.jp/get/%s"
    
    try:
        # Constructs a querystring for the RESTful URL for
        # retrieving HTML data from the relevant KEGG database entry.
        search = urlopen(URL % search_term)
        search = io.TextIOWrapper(search, encoding="UTF-8").read()
        return search
    
    except HTTPError:
        # Tries again after 10s if access is blocked.
        # This prevents the program from overloading the RESTful API.
        time.sleep(10)
        try:
            search = urlopen(URL % search_term)
            search = io.TextIOWrapper(search, encoding="UTF-8").read()
            return search

        except HTTPError:
            return None


def identify_reactions(compound):
    
    """
    Extracts the KEGG IDs of all reactions that a compound is 
    involved in from its KEGG COMPOUND database entry.
    """

    # Removes compound coefficients.
    if " " in compound:
        compound = compound.split(" ")[1]

    # Retrieves HTML data from the compound's KEGG database entry.
    search_term = "cpd:" + str(compound)
    search = get_data(search_term)
    if search is None:
        return []
    else:
        try:
            # Retrieved HTML data is read line by line.
            current_section = None
            for line in search.rstrip().split("\n"):
                section = line[:12].strip()
                if not section == "":
                    current_section = section
                    # Finds section where reactions are listed
                    # and extracts their IDs from each line.
                    if current_section == "REACTION":
                        index = search.rstrip().split("\n").index(line)
                        reactions = line[12:].split(" ")
                        for line in search.rstrip().split("\n")[int(index)+1:]:
                            # If the current line is not part of the
                            # REACTION section then the processing ends.
                            if line[:12].strip() != (""):
                                break
                            reactions_ = line[12:].split(" ")
                            reactions.extend(reactions_)

            # Eliminates erroneous IDs resulting from whitespaces.
            reactions = ['rn:' + reaction for reaction in reactions if reaction != ""]
            
            return reactions

        except UnboundLocalError:
            return []


def investigate_reaction(reaction):
    
    """
    Identifies the EC numbers of the enzymes that catalyse a reaction, along with
    the KEGG IDs of the reactants and products, from its KEGG REACTION database entry.
    """

    # Retrieves HTML data from the reaction's KEGG database entry.
    search = get_data(reaction)
    if search is None:
        return (None,)*3
    else:
        try:
            # Retrieved HTML data is read line by line.
            current_section = None
            for line in search.rstrip().split("\n"):
                section = line[:12].strip()
                if not section == "":
                    current_section = section

                    # Identifies the reaction equation and
                    # extracts the reactants and products.
                    if current_section == "EQUATION":
                        equation = line[12:]
                        components = equation.split(" <=> ") 
                        reactants_temp = components[0].split(" + ")
                        reactants = []
                        # Removes coefficients from reactants.
                        for reactant in reactants_temp:
                            if " " in reactant:
                                reactant = reactant.split(" ")[1]
                            reactants.append(reactant)
                        products_temp = components[1].split(" + ")
                        products = []
                        # Removes coefficients from products.
                        for product in products_temp:
                            if " " in product:
                                product = product.split(" ")[1]
                            products.append(product)

                    # Identifies enzymes that catalyse the
                    # reaction and extracts their EC numbers.
                    if current_section == "ENZYME":
                        if " " in line[12:]:
                            enzymes = line[12:].split(" ")
                            enzymes = ['EC:' + enzyme for enzyme in enzymes if enzyme != ""]
                        else:
                            enzymes = ["EC:" + line[12:]]

            return enzymes, reactants, products

        except UnboundLocalError:
            return (None,)*3


def chain_former(reactions, inducer, max_chain_length):
    
    """Generates linear chains of enzymes that sequentially catabolize an inducer compound."""

    products_info = {}
    reactions_info = {}
    all_chains = []

    # The IDs of unnecessary byproducts, such as H20 and NADH.
    excluded_compounds = [
    "C00035", "C00040", "C00025", "C00044", "C00007", 
    "C00005", "C00006", "C00045", "C00001", "C00010", 
    "C00024", "C00002", "C00003", "C00008", "C00009", 
    "C00011", "C00012", "C00013", "C00014", "C00019"
    ]


    def link_reactions(reaction, compound, recursion_depth=0, prior_chains=None, prior_enzymes=None, limit=max_chain_length):

        """
        Uses recursion to identify whether a reaction catabolizes a product of a reaction that 
        preceded it. Enzymes that catalyse these reactions are individually linked to generate 
        linear enzymatic chains. This process continues until chains meet the maximum chain length.
        """
    
        # Retrieves and/or stores the details of a 
        # reaction that a compound is involved in.
        if reaction not in reactions_info:
            enzymes, reactants, products = investigate_reaction(reaction)
            reaction_info = [enzymes, reactants, products]
            reactions_info[reaction] = reaction_info
        else:
            reaction_info = reactions_info[reaction]
            enzymes = reaction_info[0]
            reactants = reaction_info[1]
            products = reaction_info[2]

        if None not in reaction_info:
            starting_compound = compound
            # Reactions that catabolize the compound are processed.
            if starting_compound in reactants:
                    
                # Creates enzymatic chains by linking enzymes that catalyse the
                # reaction at the current depth to enzymes from the previous depth. 
                if (recursion_depth == 1) and (prior_enzymes is not None):
                    chains_ = [[enzyme_1, enzyme_2]
                                for enzyme_1 in prior_enzymes
                                if '-' not in enzyme_1
                                for enzyme_2 in enzymes
                                if '-' not in enzyme_2] 
                
                # Extends enzymatic chains from the previous depth.
                elif (recursion_depth > 1) and (prior_chains is not None):
                    extended_chains = [chain + [e] for e in enzymes for chain in prior_chains if '-' not in e]
                    chains_ = extended_chains
                    
                else:
                    chains_ = None
                
                # Chains at the current depth are output to terminal.
                if chains_ is not None:
                    for chain_ in chains_:
                        all_chains.append(chain_)
                        print(f"Chain identified: {' => '.join(e for e in chain_)} ")
                
                # Retrieves and/or stores the subsequent reactions 
                # that each product is involved in.
                for product in products:
                    if product not in excluded_compounds:
                        if product not in products_info:
                            reactions_2 = identify_reactions(product)
                            products_info[product] = reactions_2
                        else:
                            reactions_2 = products_info[product]

                        # Recursively either forms or extends chains based upon recursion depth.
                        if len(reactions_2) < 100:
                            for reaction_2 in reactions_2:
                                depth = recursion_depth
                                depth +=1
                                if depth == 1:
                                    link_reactions(reaction_2, product, recursion_depth=depth, prior_enzymes=enzymes)
                                elif (depth > 1) and (depth < limit):
                                    link_reactions(reaction_2, product, recursion_depth=depth, prior_chains=chains_)

    for reaction in reactions:
        link_reactions(reaction, inducer, limit=max_chain_length)

    return all_chains


def form_chains(inducer, max_chain_length):

    """
    Conducts and optimizes the chain identification procedure 
    by using multiprocessing if appropriate.
    """

    initial_reactions = identify_reactions(inducer)
    num_initial_nodes = len(initial_reactions)
    cores = multiprocessing.cpu_count()

    # Determines whether there are enough initial
    # reactions to split into multiple processes.
    if num_initial_nodes >= cores:
        processes = cores
    elif num_initial_nodes >= 2:
        processes = 2
    elif num_initial_nodes == 1:
        processes = 1
    else:
        processes = None

    if processes is not None:

        # Chains are generated in parallel depending
        # on the total number of initial reactions.
        if processes >= 2:
            initial_reactions = np.array(initial_reactions, dtype=object)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                futures = []
                for data in np.array_split(initial_reactions, processes):
                    future = executor.submit(chain_former, data, inducer, max_chain_length)
                    futures.append(future)

                chains_all = []
                for future in futures:
                    result = future.result()
                    for chain in result:
                        if chain not in chains_all:
                            chains_all.append(chain)

        # Chains are generated sequentially if
        # there is only one initial reaction.
        else:
            chains_all = []
            chains_unfiltered = chain_former(initial_reactions, inducer, max_chain_length)
            for chain in chains_unfiltered:
                if chain not in chains_all:
                    chains_all.append(chain)
    else:
        chains_all = []

    if len(chains_all) > 0:
        return chains_all
    else:
        sys.exit(f"No chains were identified for '{inducer}'")


def retrieve_genes(chain):
    
    """Identifies organisms which possess all enzymes within a chain and retrieves the corresponding genes."""

    all_genes = {}

    for n in range(len(chain)):

        # Data is retrieved from the enzyme's KEGG database entry.
        enzyme = chain[n].lower()
        search = get_data(enzyme)

        # Data is read line by line
        if search is not None:
            x = []
            current_section = None
            for line in search.rstrip().split("\n"):
                section = line[:12].strip()
                if not section == "":
                    current_section = section

                    # Identifies the organisms that possess the enzyme
                    # and the genes that encode it in their genomes.
                    if current_section == "GENES":
                        index = search.rstrip().split("\n").index(line)
                        x.append(line[12:].split(": "))
                        for line in search.rstrip().split("\n")[int(index)+1:]:
                            if line[:12].strip() != "":
                                break
                            x.append(line[12:].split(": "))

            # Organisms and their genes are stored in a dataframe.
            all_genes["d" + str(n)] = pd.DataFrame(x, columns=["Organism", str(enzyme) + "_gene(s)"])
    
    gene_sets = len(all_genes)
    if gene_sets > 1:    
        try:
            # Merges the dataframes to form one that only includes 
            # organisms that possess all enzymes within the chain.
            for x in range(gene_sets-1):
                all_genes["d" + str(x+1)] = pd.merge(left = all_genes["d" + str(x)], right = all_genes["d" + str(int(x)+1)], left_on = "Organism", right_on = "Organism")

            filtered_genes = all_genes["d" + str(max(range(gene_sets)))] 

            return filtered_genes

        except KeyError:
            return None


class Biosensor(typ.NamedTuple):
    
    """A constructor for storing predicted biosensors."""

    operon: str
    regulator: str
    regulator_score: int
    regulator_annotation: str
    organism_code: str
    genes: dict
    gene_positions: dict


def predict_biosensors(df, genome_assemblies, genome_files):
    
    """
    Applies the in_regulon function to each row of a dataframe of organisms
    and their genes that encode each enzyme within a chain. 
    """
    
    biosensors = []

    def in_regulon(row, columns):
        
        """
        Determines whether an organism's genes that encode enzymes within a chain have 
        genetic organisations that are characteristic of operons and applies the identify_regulator 
        function to each predicted operon to identify their putative transcriptional regulators.
        """
        
        # Uses the organism code of a gene set to
        # identify the correct genome assembly to use.
        organism_code = row[0]
        index = genome_assemblies.index[genome_assemblies["Organism code"] == str(organism_code).lower()]
        if not index.empty:
            assembly = genome_assemblies["Assembly"][index].to_string(index=False, header=False)
            
            # Finds the feature table genome that that has 
            # the relevant genome assembly in its filename.
            match = [s for s in genome_files if assembly[3:] in s]
        
            if match != []:

                # Reads and parses the correct feature table genome.
                genome = pd.read_csv(match[0], sep="\t")
                genome.drop(genome[genome["# feature"] == "gene"].index, inplace=True)
                genome = genome.reset_index(drop=True)
                num_cols = len(columns)

                gene_positions = {}
                genes_ = {}
                for n in range(1, num_cols):
                    genes = row[n]
                    genes_[n] = genes

                # Genes that encode the first enzyme within a chain
                # are used as starting genes for potential operons.
                starting_genes = genes_[1]
                starting_genes = starting_genes.split(" ")
                all_genes = [row[n].split(" ") for n in range(1, num_cols)]
                all_genes = list(chain(*all_genes))
                avoid_repeats = []

                for x in range(len(starting_genes)):
                    starting_gene = starting_genes[x]
                    starting_gene_ = starting_gene.split("(")[0]
                    operon = [starting_gene]
                    avoid_repeats.append(x)

                    try:
                        # Determines the index position and strand orientation of the starting gene.
                        starting_position = genome.index[genome["locus_tag"] == starting_gene_].tolist()[0]
                        starting_orientation = genome["strand"][starting_position]

                        # Determines the index position and strand orientation of all other genes.
                        for y in range(len(all_genes)):
                            if y not in avoid_repeats:
                                gene = all_genes[y]
                                gene_ = gene.split("(")[0]
                                position = genome.index[genome["locus_tag"] == gene_].tolist()[0]
                                gene_orientation = genome["strand"][position]

                                try:
                                    # Genes that are nearby the starting gene and have 
                                    # the same strand orientation are placed in the operon.
                                    if (abs(int(starting_position) - int(position)) < 30) and (starting_orientation == gene_orientation):
                                        operon.append(gene)
                                        if y <= len(starting_genes)-1:
                                            avoid_repeats.append(y)
                                        if gene not in gene_positions:
                                            gene_positions[gene] = position

                                except ValueError:
                                    pass
                        
                        # Identifies putative regulators of predicted operons
                        # to predict potential biosensors for the inducer compound.
                        if len(operon) > 1:
                            if starting_gene not in gene_positions:
                                gene_positions[starting_gene] = starting_position
                            regulator, score, annotation = identify_regulator(genome, operon, starting_orientation, gene_positions)
                            if None not in [regulator, score, annotation]:
                                biosensor = Biosensor(operon, regulator, score, annotation, organism_code, genes_, gene_positions)
                                biosensors.append(biosensor)
                    
                    except IndexError:
                        pass

    df.apply(lambda x: in_regulon(x, df.columns), axis=1)

    return biosensors


def identify_regulator(genome, operon, operon_orientation, gene_positions):
    
    """
    Identifies putative regulators of predicted operons based upon a conceptual model.
    Regulatory genes are often situated directly upstream of their corresponding operon
    and on the opposite DNA strand.
    """

    positions = [gene_positions[gene] for gene in operon]
    
    try:
        # Selects regulators situated on the reverse DNA strand
        # that are upstream of an operon on the forward DNA strand.
        if operon_orientation == "+":
            start_position = min(positions)
            start_seqtype = genome["seq_type"][min(positions)]
            regulators = genome[genome["name"].str.contains("regulator|repressor|activator") == True]
            regulators = regulators[regulators["strand"] == "-"]
            regulators = regulators[regulators["seq_type"] == start_seqtype]
            reg_positions = regulators.index.to_list()
            reg_positions = [reg_position for reg_position in reg_positions if reg_position <= start_position]

        # Selects regulators situated on the forward DNA strand
        # that are upstream of an operon on the reverse DNA strand.
        elif operon_orientation == "-":
            start_position = max(positions)
            start_seqtype = genome["seq_type"][max(positions)]
            regulators = genome[genome["name"].str.contains("regulator|repressor|activator") == True]
            regulators = regulators[regulators["strand"] == "+"]
            regulators = regulators[regulators["seq_type"] == start_seqtype]
            reg_positions = regulators.index.to_list()
            reg_positions = [reg_position for reg_position in reg_positions if reg_position >= start_position]

        else:
            reg_positions = None
    
    except (AttributeError, IndexError):
        reg_positions = None

    if (reg_positions is not None) and (len(reg_positions) > 0):
        try:
            # Finds the shortest distance between the operon
            # and the previously selected regulators.
            arr = np.asarray(reg_positions)
            distances = abs(arr - start_position)
            min_distance = np.amin(distances)
            idx = int(np.where(distances == min_distance)[0][0])

            # Finds the closest of these regulators to the operon.
            regulator_position = reg_positions[idx]
            closest_regulator = regulators["locus_tag"][regulator_position]
            annotation =  regulators["name"][regulator_position]
            
            score = 0

            # The closest regulator does not have any points deducted
            # from its score if it situated right next to the operon.
            if abs(regulator_position - start_position) == 1:
                pass

            # Scores closest regulators that are situated on the forward
            # DNA strand and do not directly neighbour the operon.
            elif regulator_position > start_position:
                for n in range(1, min_distance):
                    orient = genome["strand"][start_position + n]
                    if orient != operon_orientation:
                        score -= 2
                    else:
                        score -= 1

            # Scores closest regulators that are situated on the reverse
            # DNA strand and do not directly neighbour the operon.
            elif regulator_position < start_position:
                for n in range(1, min_distance):
                    orient = genome["strand"][regulator_position + n]
                    if orient != operon_orientation:
                        score -= 2
                    else:
                        score -= 1
                
            else:
                score = None
            
            return closest_regulator, score, annotation

        except (ValueError, TypeError):
            return (None,)*3
    else:
        return (None,)*3
    

def process_chain(chain, inducer, genome_assemblies, genome_files):

    """
    Applies the biosensor prediction algorithms to an enzymatic chain 
    and outputs the data to csv files in an understandable format.
    """

    # Identifies whether there are any genomes that
    # encode all enzymes within the chain.
    df = retrieve_genes(chain)    
    if df is not None:
        if (not df.empty) and (len(df.columns) > 2):

            # Determines whether multiprocessing should be used on the data
            # and calculates the number of processes to conduct.
            cores = multiprocessing.cpu_count()    
            if len(df) >= cores:
                processes = cores
            elif len(df) >= 2:
                processes = 2
            else:
                processes = None
            
            # Splits the data evenly and applies the biosensor prediction 
            # algorithms to the data using multiprocessing.
            if processes is not None:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    futures = [executor.submit(predict_biosensors, data, genome_assemblies, genome_files) for data in np.array_split(df, processes)]            
                    results = [future.result() for future in futures]
                    regulons = [inner for outer in results for inner in outer]

            # If there is not enough data for multiprocessing,
            # the data is instead processed sequentially.
            else:   
                regulons = predict_biosensors(df, genome_assemblies, genome_files)
            
            # If biosensors were predicted, they are ranked in order
            # of their scores and formatted for data output.
            num_regulons = len(regulons)
            if num_regulons > 0:
                regulons.sort(key=lambda x: x.regulator_score, reverse=True)
                chain_length = len(chain)

                # Prepares output directory names for the data.
                root = str(inducer) + "_results"
                subroot = "chainlength=" + str(chain_length)
                cols = df.columns.values.tolist()
                name = inducer + "_" + "_".join(str(col)[3:9] for col in cols[1:]) + ".csv"
                
                # Formats the data for csv files.
                header = ["Organism_code"] + [cols[x] for x in range(1, chain_length+1)] + ["Operon", 
                            "Regulator", 
                            "Regulator_score", 
                            "Regulator_annotation"]
                
                data = [[regulon.organism_code] + 
                        [regulon.genes[x] for x in range(1, chain_length+1)] +
                        [" ".join(str(gene) for gene in regulon.operon),
                        regulon.regulator,
                        regulon.regulator_score,
                        regulon.regulator_annotation] 
                        for regulon in regulons]
                
                # Creates directories and outputs
                # data as .csv files to them.
                try:
                    path = os.path.join(root, subroot, name)
                    with open(path, "w", encoding="UTF8", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow(header)
                        writer.writerows(data)

                except FileNotFoundError:
                    dir = os.path.join(root, subroot)
                    os.makedirs(dir)
                    path = os.path.join(root, subroot, name)

                    with open(path, "w", encoding="UTF8", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow(header)
                        writer.writerows(data)

                return num_regulons


def main():
    
    t1 = time.time()
    args = argument_parser()
    inducer = args.compound
    max_chain_length = args.length

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
    
    print("{}Identifying enzymatic chains for '{}' with maximum chain length set to {}...{}".format("\n", inducer, max_chain_length, "\n"))
    chains = form_chains(inducer, max_chain_length)
    total_chains = len(chains)
    print("{}{} unique chains were identified.".format("\n", total_chains))
    pd.set_option('mode.chained_assignment', None)
    print("{}Processing {} chains...{}".format("\n", total_chains, "\n"))
    
    total_biosensors = 0
    for n in tqdm(range(total_chains)):
        num_regulons = process_chain(chains[n], inducer, genome_assemblies, genome_files)
        if num_regulons is not None:
            total_biosensors += num_regulons
    t2 = time.time()
    if total_biosensors > 0:
        print("{}Processing is complete. {} potential biosensors were identified for '{}'. Results have been deposited to '{}_results'. Total runtime: {}s.".format("\n", total_biosensors, inducer, inducer, round(t2-t1, 2)))
    else:
        print("{}Processing is complete. {} potential biosensors were identified for '{}'. Total runtime: {}".format("\n", total_biosensors, inducer, round(t2-t1, 2)))


if __name__ == "__main__":
    main()
