"""
Provides a command-line interface.
"""

import argparse


def argument_parser():
    """
    Parses command-line arguments and provides usage information.
    """
    parser = argparse.ArgumentParser(
        prog="TFBMiner",
        usage="py -m TFBMiner [-h] [-l L] [-s S] compound",
        description = "TFBMiner: Identifies putative transcription factor-based biosensors for a given compound."
    )
    parser.add_argument(
        "compound", 
        type=str, 
        help="Enter the KEGG Compound ID of the inducer compound."
    )
    parser.add_argument(
        "-l", 
        "--max_chain_length", 
        type=int, 
        help="Enter the maximum length for identified enzymatic chains.", 
        default=3
    )
    parser.add_argument(
        "-s", 
        "--single_gene_operons", 
        help="Choose whether to predict biosensors for rare, potential single-gene operons (y/n).", 
        default="n"
    )
    parser.add_argument(
        "-g",
        "--genome_files_path",
        help="Enter the absolute path of the 'genome_files' directory that contains feature table genomes of bacteria held on the KEGG GENOME database. If unspecified, TFBMiner will try to locate it within the user's home directory.",
        default="unspecified"
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
        help="Enter the absolute path of the desired output directory. If unspecified, TFBMiner will output the results to the user's home directory.",
        default="unspecified"
    )

    args = parser.parse_args()
    return args