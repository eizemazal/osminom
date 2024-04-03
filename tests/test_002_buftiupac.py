import pytest

from buftinom.iupac import Iupac
from buftinom.smiles_parser import SmilesParser


@pytest.fixture
def parser():
    return SmilesParser(debug=True)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C", "methane"),
        ("CC", "ethane"),
        ("C=C", "ethene"),
        ("C#C", "ethyne"),
        #
        ("CCC", "propane"),
        ("C=CC", "propene"),
        ("C#CC", "propyne"),
        #
        ("CCCC", "butane"),
        ("C=CCC", "butene"),
        ("CCC=C", "butene"),
        ("C#CCC", "butyne"),
        ("CCC#C", "butyne"),
        ("CC=CC", "but-2-ene"),
        ("CC#CC", "but-2-yne"),
        #
        ("CCCCC", "pentane"),
        ("CC=CCC#CC", "hept-2-en-5-yne"),
        ("CCC=CC#CC", "hept-4-en-2-yne"),
        #
        ("C=C#C", "propen-2-yne"),
    ],
)
def test(parser, smiles, expected):
    (mol,) = parser.parse(smiles)
    iupac = Iupac(mol)

    assert iupac.simple_chain_name(mol._atoms) == expected
