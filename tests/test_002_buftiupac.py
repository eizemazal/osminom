import pytest

from buftinom.iupac import Iupac, iupac2str
from buftinom.smileg import debug_atoms
from buftinom.smiles_parser import SmilesParser

debug_atoms(True)


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
        ("CC=C=CC", "pent-2,3-diene"),
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
        ("CCCCC(C(C)C)CCCCC", "5-(1-methylethyl)decane"),
        ("CCC(C)C(C(C)CC)CCC", "3,5-dimethyl-4-propylheptane"),
        ("CCC(C)C(C(C)C)C=C", "3-(1-methylethyl)-4-methylhexene"),
        ("CCC(C)C(C(C)CC)C=C", "4-ethenyl-3,5-dimethylheptane"),
        ("CCC(C)C(C)CCCCC(C)C", "2,7,8-trimethyldecane"),
        ("CCCC(C)C(CC)CCC", "4-ethyl-5-methyloctane"),
        ("CC(C)CCC(CC)CC", "5-ethyl-2-methylheptane"),
        ("CCC=CCC(C)C", "6-methylhept-3-ene"),
        ("CC=CCC#CC", "hept-2-en-5-yne"),
        ("C(C)(C)C(C)(C)", "2,3-dimethylbutane"),
        ("C(C)(C)C(C)(C)CCC", "2,3,3-trimethylhexane"),
        ("CCCCC(C(C)C(C)C)(C(C)C(C)C)CCCC", "5,5-bis(1,2-dimethylpropyl)nonane"),
        ("CCCC(C(C)(C)C)CC", "3-ethyl-2,2-dimethylhexane"),
        ("CCCCC(C(C)C(C)CC)C(CC)CCC", "5-butyl-6-ethyl-3,4-dimethylnonane"),
        ("CCCCCC(C(C)C(C)CC)C(CC)CCCC", "6-(1,2-dimethylbutyl)-5-ethylundecane"),
    ],
)
def test_various_molecules(parser, smiles, expected):
    (mol,) = parser.parse(smiles)
    iupac = Iupac(mol)

    mol.print_table()
    iupac.decomposition.print()

    name = iupac.decompose_name(iupac.decomposition)
    assert iupac2str(name) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC1CCCCC1", "1-methylcyclohexane"),
        ("C1CC(C)CCC1", "1-methylcyclohexane"),
        ("C1CCCCC1C", "1-methylcyclohexane"),
        ("CC1C(C)CC(C)CC1", "1,2,4-trimethylcyclohexane"),
    ],
)
def test_cycle_names(parser, smiles, expected):
    (mol,) = parser.parse(smiles)
    iupac = Iupac(mol)

    mol.print_table()
    iupac.decomposition.print()

    name = iupac.decompose_name(iupac.decomposition)
    assert iupac2str(name) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCO", "propanol"),
        ("OCCC", "propanol"),
        ("CC(O)C", "propan-2-ol"),
        ("CC(O)C(O)C", "butan-2,3-diol"),
        ("CCC(=O)O", "propanoic acid"),
        ("CCC(O)=O", "propanoic acid"),
        ("CCCN", "propanamine"),
        ("CC(C)CCCO", "4-methylpentanol"),
        ("C1CCC(O)CC1", "cyclohexanol"),
    ],
)
def test_functional_group_naming(parser, smiles, expected):
    (mol,) = parser.parse(smiles)
    iupac = Iupac(mol)

    mol.print_table()
    iupac.decomposition.print()
    for f in iupac.alg.functional_groups.items():
        print(f)

    name = iupac2str(iupac.decompose_name(iupac.decomposition))
    print(name)
    assert name == expected
