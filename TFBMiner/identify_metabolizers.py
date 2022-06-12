"""
Conducts the chain identification stage of the pipeline.
"""

from itertools import chain
import concurrent.futures
import multiprocessing
import sys

import numpy as np

from TFBMiner import acquire_data


class MetabolizerIdentifier:
    
    def __init__(self, inducer):
        self.inducer = inducer

    def identify_chains(self, reactions, max_chain_length):
        """
        Generates linear chains of enzymes that sequentially catabolize an inducer compound.
        """
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

        def link_reactions(reaction, compound, recursion_depth=0, limit=max_chain_length, prior_chains=None, prior_enzymes=None):
            """
            Uses recursion to identify whether a reaction catabolizes a product of a reaction that 
            preceded it. Enzymes that catalyse these reactions are individually linked to generate 
            linear enzymatic chains. This process continues until chains meet the maximum chain length.
            """
            # Retrieves and/or stores the details of a 
            # reaction that a compound is involved in.
            if reaction not in reactions_info:
                enzymes, reactants, products = acquire_data.reaction_details(reaction)
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
                                reactions_2 = acquire_data.identify_reactions(product)
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
            link_reactions(reaction, self.inducer)

        return all_chains


    def execute_chain_identification(self, max_chain_length):
        """
        Conducts and optimizes the chain identification procedure 
        by using multiprocessing if appropriate.
        """
        print("{}Identifying enzymatic chains for {} with maximum chain length set to {}...{}".format("\n", self.inducer, max_chain_length, "\n"))
        reactions = acquire_data.identify_reactions(self.inducer)
        total_reactions = len(reactions)
        cores = multiprocessing.cpu_count()

        processes = cores if total_reactions>=cores else 2 if total_reactions>=2 else 1 if total_reactions==1 else None
        if processes is not None:
            if processes >= 2:
                # Identifies enzymatic chains concurrently.
                reactions = np.array(reactions, dtype=object)
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    futures = []
                    for data in np.array_split(reactions, processes):
                        future = executor.submit(self.identify_chains, data, max_chain_length)
                        futures.append(future)
                    
                    chains = []
                    for future in futures:
                        result = future.result()
                        for chain in result:
                            if chain not in chains:
                                chains.append(chain)
            else:
                chains = []
                chains_unfiltered = self.identify_chains(reactions, max_chain_length)
                for chain in chains_unfiltered:
                    if chain not in chains:
                        chains.append(chain)
        else:
            chains = []
        
        total_chains = len(chains)
        if total_chains > 0:
            print("{}{} unique chains were identified.".format("\n", total_chains))
            return chains, total_chains
        else:
            sys.exit(f"No chains were identified for {self.inducer}")


    def identify_single_metabolizers(self, reactions):
        """
        Finds only initial enzymes that metabolize a compound.
        """
        enzymes = []
        for reaction in reactions:
            enzymes_, reactants, products = acquire_data.reaction_details(reaction)
            if self.inducer in reactants:
                metabolizers = [enzyme_ for enzyme_ in enzymes_ if "-" not in enzyme_]
                if len(metabolizers) > 0:
                    enzymes += metabolizers
                    for enzyme in metabolizers:
                        print(f"Metabolizer identified: {enzyme}")
        return enzymes


    def execute_single_metabolizer_identification(self):
        """
        Conducts and optimizes the single enzyme metabolizer 
        identification procedure by using multiprocessing if appropriate.
        """
        print("{}Identifying single enzymes that metabolize {}...{}".format('\n', self.inducer, '\n'))
        reactions = acquire_data.identify_reactions(self.inducer)
        total_reactions = len(reactions)
        cores = multiprocessing.cpu_count()
        
        processes = cores if total_reactions>=cores else 2 if total_reactions>=2 else 1 if total_reactions==1 else None
        if processes is not None:
            if processes >= 2:
                # Identifies single enzyme metabolizers concurrently.
                reactions = np.array(reactions, dtype=object)
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    futures = [executor.submit(self.identify_single_metabolizers, data) for data in np.array_split(reactions, processes)]            
                    enzymes = [future.result() for future in futures]
            elif processes == 1:
                enzymes = self.identify_single_metabolizers(reactions[0])
        else:
            enzymes = []

        # Removes repeated enzymes.
        enzymes = list(chain(*enzymes))
        enzymes = list(dict.fromkeys(enzymes))
        total_enzymes = len(enzymes)
        if total_enzymes > 0:
            print("{}{} unique enzymes were identified as metabolizers of {}.".format('\n', total_enzymes, self.inducer))
            return enzymes, total_enzymes
        else:
            sys.exit(f"No chains were identified for {self.inducer}")
