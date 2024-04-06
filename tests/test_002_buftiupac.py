import pytest

from buftinom.iupac import Iupac, iupac2str
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
def test_simple_chain_name(parser, smiles, expected):
    (mol,) = parser.parse(smiles)
    iupac = Iupac(mol)
    name = iupac.decompose_name(iupac.decomposition)

    assert iupac2str(name) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC(C)(C)CC", "3,3-dimethylpentane"),
        ("CC(C)CC(C)C", "2,4-dimethylpentane"),
        # Note here he found that it is hexane, ill provide numeration for clearance
        #  3 21 4   56
        ("CC(CC)C(C)CC", "3,4-dimethylhexane"),
        ("CCCCC(CC)C(C)CCC", "5-ethyl-4-methylnonane"),
    ],
)
def test_decomposed_multichain(parser, smiles, expected):
    (mol,) = parser.parse(smiles)
    iupac = Iupac(mol)

    mol.print_table()
    iupac.decomposition.print()

    name = iupac.decompose_name(iupac.decomposition)
    assert iupac2str(name) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC(C)(C)CCC", "3,3-dimethylhexane"),
        ("CCC(C)CC(C)CC", "3,5-dimethylheptane"),
        ("CCC(C)C(C(C)CC)CCC", "3,5-dimethyl-4-propylheptane"),
        ("CCC(C)C(C(C)C)C=C", "3-(1-methylethyl)-4-methylhexene"),
        ("CCC(C)C(C(C)CC)C=C", "4-ethenyl-3,5-dimethylheptane"),
        ("CCC(C)C(C)CCCCC(C)C", "2,7,8-trimethyldecane"),
        ("CCCC(C)C(CC)CCC", "4-ethyl-5-methyloctane"),
    ],
)
def test_various_molecules(parser, smiles, expected):
    (mol,) = parser.parse(smiles)
    iupac = Iupac(mol)

    mol.print_table()
    iupac.decomposition.print()

    name = iupac.decompose_name(iupac.decomposition)
    assert iupac2str(name) == expected
