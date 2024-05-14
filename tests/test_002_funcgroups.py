import pytest

from osminom.structure2name.funcgroups import acid_matcher, alco_matcher
from osminom.structure2name.molecule import Molecule


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
    mol = Molecule.from_smiles(smiles)
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
    mol = Molecule.from_smiles(smiles)
    acid = acid_matcher(mol)

    nmatches = 0
    for atom in mol.atoms:
        mtch = acid.matches(atom)
        if mtch is None:
            continue

        print(mtch)
        nmatches += 1

    assert nmatches == expected_nmatches
