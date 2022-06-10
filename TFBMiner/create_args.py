"""
Instantiates command-line arguments to control TFBMiner.
"""

import argparse


def argument_parser():
    """
    Creates command-line options for the user.
    """
    parser = argparse.ArgumentParser(
        prog="TFBMiner",
        usage="py -m TFBMiner [-h] [-l L] [-s S] compound",
        description = "TFBMiner: Identifies putative transcription factor-based biosensors for a given compound.")
    parser.add_argument(
        "compound", 
        type=str, 
        help="Enter the KEGG Compound ID of the inducer compound.")
    parser.add_argument(
        "-l", 
        "--length", 
        type=int, 
        help="Enter the maximum length of the enzymatic chains.", 
        default=3)
    parser.add_argument(
        "-s", 
        "--single_gene_operons", 
        help="Choose whether to predict biosensors for rare, potential single-gene operons (y/n).", 
        default="n")

    args = parser.parse_args()
    return args