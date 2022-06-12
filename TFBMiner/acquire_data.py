from urllib.error import HTTPError
from urllib.request import urlopen
import time
import io


def get_data(search_term):
    """
    Retrieves data held in KEGG database entries using KEGG's RESTful API.
    """
    URL = "http://rest.kegg.jp/get/%s"
    
    try:
        # Constructs and uses a querystring for the RESTful URL to
        # retrieve HTML data for the relevant KEGG database entry.
        data = urlopen(URL % search_term)
        data = io.TextIOWrapper(data, encoding="UTF-8").read()
        return data
    
    except HTTPError:
        # Tries again after 10s if access is blocked.
        # This prevents the program from overloading the RESTful API.
        time.sleep(10)
        try:
            data = urlopen(URL % search_term)
            data = io.TextIOWrapper(data, encoding="UTF-8").read()
            return data

        except HTTPError:
            return None


def identify_reactions(compound):
    """
    Extracts the KEGG IDs of all reactions that a compound is 
    involved in from its KEGG COMPOUND database entry.
    """
    # Removes coefficients that obfuscate compound IDs.
    if " " in compound:
        compound = compound.split(" ")[1]
    search_term = "cpd:" + str(compound)
    data = get_data(search_term)
    if data is None:
        return []
    else:
        try:
            # Retrieved HTML data is read line by line.
            current_section = None
            for line in data.rstrip().split("\n"):
                section = line[:12].strip()
                if not section == "":
                    current_section = section
                    # Finds section where reactions are listed
                    # and extracts their IDs from each line.
                    if current_section == "REACTION":
                        index = data.rstrip().split("\n").index(line)
                        reactions = line[12:].split(" ")
                        for line in data.rstrip().split("\n")[int(index)+1:]:
                            # If the current line is not part of the REACTION 
                            # section then the processing ends.
                            if line[:12].strip() != (""):
                                break
                            reactions_ = line[12:].split(" ")
                            reactions.extend(reactions_)
            # Eliminates erroneous IDs resulting from whitespaces.
            reactions = ['rn:' + reaction for reaction in reactions if reaction != ""]
            
            return reactions

        except UnboundLocalError:
            return []


def reaction_details(reaction):
    """
    Identifies EC numbers of the enzymes that catalyse a reaction, along with
    the KEGG IDs of the reactants and products, from its KEGG database entry.
    """
    data = get_data(reaction)
    if data is None:
        return (None,)*3
    else:
        try:
            # Retrieved HTML data is read line by line.
            current_section = None
            for line in data.rstrip().split("\n"):
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
            

def retrieve_encoders(enzyme):
    """
    Finds the genes and organisms that encode an enzyme.
    """
    data = get_data(enzyme)
    if data is not None:
        encoders = []
        current_section = None
        for line in data.rstrip().split("\n"):
            section = line[:12].strip()
            if not section == "":
                current_section = section
                if current_section == "GENES":
                    index = data.rstrip().split("\n").index(line)
                    encoders.append(line[12:].split(": "))
                    for line in data.rstrip().split("\n")[int(index)+1:]:
                        if line[:12].strip() != "":
                            break
                        encoders.append(line[12:].split(": "))
        return encoders