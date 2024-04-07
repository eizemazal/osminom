import pytest

from buftinom.funcgroups import acid_matcher2, alco_matcher
from buftinom.smiles_parser import SmilesParser


@pytest.mark.parametrize(
    "smiles, expected_nmatches",
    [
        ("C", 0),
        ("O", 0),
        ("CCC", 0),
        ("OCCC", 1),
        ("CC(O)C", 1),
        ("CCCO", 1),
        ("CC(O)CCO", 2),
    ],
)
def test_alcohol(smiles, expected_nmatches):
    (mol,) = SmilesParser().parse(smiles)
    alcohol = alco_matcher(mol)

    nmatches = 0
    for atom in mol.atoms:
        mtch = alcohol.matches(atom)
        if mtch is None:
            continue

        print(mtch)
        nmatches += 1

    assert nmatches == expected_nmatches


@pytest.mark.parametrize(
    "smiles, expected_nmatches",
    [
        ("C", 0),
        ("O", 0),
        ("CCC", 0),
        ("C(=O)O", 1),
        ("CC(=O)O", 1),
        ("CCCC(=O)O", 1),
    ],
)
def test_acid(smiles, expected_nmatches):
    (mol,) = SmilesParser().parse(smiles)
    acid = acid_matcher2(mol)

    nmatches = 0
    for atom in mol.atoms:
        mtch = acid.matches(atom)
        if mtch is None:
            continue

        print(mtch)
        nmatches += 1

    assert nmatches == expected_nmatches
