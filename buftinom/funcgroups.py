from buftinom.smiles_parser import SmilesParser


def alcohol():
    (mol,) = SmilesParser().parse("CO")
    return mol


def carboxylic_acid():
    (mol,) = SmilesParser().parse("C(=O)O")
    return mol
