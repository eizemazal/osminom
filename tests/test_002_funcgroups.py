from buftinom.funcgroups import alcohol, carboxylic_acid
from buftinom.smiles_parser import SmilesParser


def test_alcohol():
    alc = alcohol()
    alc.print_table()

    (mol,) = SmilesParser().parse("C(O)C")

    oxy = mol.atom("O")
    mol.neighbors

    mol.print_table()


def test_acid():
    acid = carboxylic_acid()
    acid.print_table()

    (mol,) = SmilesParser().parse("C(=O)(O)CCC")
    mol.print_table()
