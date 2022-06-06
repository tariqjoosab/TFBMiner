"""
Contains a class constructor for storing predicted biosensors.
"""

import typing as typ


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