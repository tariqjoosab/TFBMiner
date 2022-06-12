"""
Contains functions for predicting operons and biosensors within specific genomes.
"""

from itertools import chain
import multiprocessing
import concurrent.futures
import typing as typ

import numpy as np
import pandas as pd


class Biosensor(typ.NamedTuple):
    """
    Stores biosensors and their relevant attributes in named tuples.
    """
    operon: str
    regulator: str
    regulator_score: int
    regulator_annotation: str
    organism_code: str
    genes: dict
    gene_positions: dict


def select_genome(organism_code, genome_assemblies, genome_files):
    """
    Identifies and returns an organism's GenBank genome assembly.
    """
    # Finds the relevant genome assembly for an organism.
    index = genome_assemblies.index[genome_assemblies["Organism code"] == str(organism_code).lower()]
    if not index.empty:
        assembly = genome_assemblies["Assembly"][index].to_string(index=False, header=False)
        
        # Finds the feature table genome that has 
        # the relevant genome assembly in its filename.
        match = [s for s in genome_files if assembly[3:] in s]
    
        if match != []:

            # Reads and parses the correct feature table genome.
            genome = pd.read_csv(match[0], sep="\t")
            genome.drop(genome[genome["# feature"] == "gene"].index, inplace=True)
            genome = genome.reset_index(drop=True)

            return genome


def identify_regulator(genome, operon, operon_orientation, gene_positions):
    """
    Identifies putative regulators of predicted operons based upon a conceptual model;
    regulatory genes are often situated directly upstream of their corresponding operon
    and on the opposite DNA strand.
    """
    positions = [gene_positions[gene] for gene in operon]

    try:
        # Selects regulators situated on the reverse DNA strand
        # that are upstream of an operon on the forward DNA strand.
        if operon_orientation == "+":
            start_position = min(positions)
            start_seqtype = genome["seq_type"][start_position]
            regulators = genome[genome["name"].str.contains("regulator|repressor|activator") == True]
            regulators = regulators[regulators["strand"] == "-"]
            regulators = regulators[regulators["seq_type"] == start_seqtype]
            reg_positions = regulators.index.to_list()
            reg_positions = [reg_position for reg_position in reg_positions if reg_position <= start_position]

        # Selects regulators situated on the forward DNA strand
        # that are upstream of an operon on the reverse DNA strand.
        elif operon_orientation == "-":
            start_position = max(positions)
            start_seqtype = genome["seq_type"][start_position]
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


def predict_biosensors(df, genome_assemblies, genome_files, single_gene_operons=False):
    """
    Applies a regulon identification algorithm to each row of a dataframe 
    wherein each row contains a complete set of genes that encode an enzymatic 
    chain or single enzyme, within a specific organism.
    """
    biosensors = []

    def identify_regulons(row, columns):
        """
        Determines whether genes encoding a chain are organized in a manner
        that is characteristic of an operon, wherein they are clustered together 
        on the same DNA strand, and, if so, identifies a potential transcriptional 
        regulator that could be a biosensor for the compound that the operon metabolizes.
        """
        organism_code = row[0]
        genome = select_genome(organism_code, genome_assemblies, genome_files)

        if genome is not None:
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

    def identify_single_gene_regulons(row):
        """
        Identifies possible single-gene operons and their likely
        transcriptional regulators to predict potential biosensors.
        """
        organism_code = row[0]
        genome = select_genome(organism_code, genome_assemblies, genome_files)
        if genome is not None:
            genes = row[1].split(" ")
            for gene in genes:
                gene = gene.split("(")[0]
                operon = [gene]
                try:
                    gene_position = genome.index[genome["locus_tag"] == gene].tolist()[0]
                    gene_positions = {gene: gene_position}
                    gene_orientation = genome["strand"][gene_position]
                    regulator, score, annotation = identify_regulator(genome, operon, gene_orientation, gene_positions)
                    if None not in [regulator, score, annotation]:
                        biosensor = Biosensor(operon, regulator, score, annotation, organism_code, {1: gene}, gene_positions)
                        biosensors.append(biosensor)
                except IndexError:
                    pass

    if single_gene_operons==False:
        df.apply(lambda x: identify_regulons(x, df.columns), axis=1)

    elif single_gene_operons==True:
        df.apply(lambda x: identify_single_gene_regulons(x), axis=1)

    return biosensors


def execute_biosensor_predictions(df, genome_assemblies, genome_files, single_gene_operons=False):
    """
    Optimizes the execution of the biosensor prediction functions
    based upon the size of the dataframe.
    """
    # Determines whether the data is large enough for multiprocessing 
    # and calculates the number of processes to conduct.
    total_organisms = len(df)
    cores = multiprocessing.cpu_count()
    if total_organisms >= cores:
        processes = cores
    elif total_organisms >= 2:
        processes = 2
    else:
        processes = None
    
    # Splits the data and conducts multiple processes.
    if processes is not None:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(predict_biosensors, data, genome_assemblies, genome_files, single_gene_operons=single_gene_operons) for data in np.array_split(df, processes)]            
            results = [future.result() for future in futures]
            biosensors = [inner for outer in results for inner in outer]
    # Otherwise sequentially processes the data  
    else: 
        biosensors = predict_biosensors(df, genome_assemblies, genome_files)

    return biosensors