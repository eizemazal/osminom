from contextlib import contextmanager

import pytest

from buftinom.smileg import Bond, BondType, MoleculeConstructor
from buftinom.smiles_lexer import SmilesLexer
from buftinom.smiles_parser import SmilesParser
from tests import test_data


@contextmanager
def maybe_raises(exc, should_raise):
    if not should_raise:
        yield
        return

    with pytest.raises(exc):
        yield


@pytest.mark.parametrize(
    "smiles,valid",
    test_data.smiles,
)
def test_lexer(smiles, valid):
    lexer = SmilesLexer()

    with maybe_raises(ValueError, not valid):
        lexer.test(smiles)


@pytest.mark.parametrize(
    "smiles,valid",
    test_data.smiles + test_data.parser_tests,
)
def test_parser(smiles, valid):
    parser = SmilesParser(debug=True)

    with maybe_raises(ValueError, not valid):
        cosmos = parser.parse(smiles)

        for mol in cosmos:
            print(mol.print_table())


@pytest.fixture
def parser():
    return SmilesParser(debug=True)


def test_many_molecules(parser):
    assert len(parser.parse("C.C")) == 2


def bond_values(mol: MoleculeConstructor):
    return [b for _, b, _ in mol.unuque_bonds()]


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C-C", [Bond(BondType.SINGLE)]),
        ("C=C", [Bond(BondType.DOUBLE)]),
        ("C#C", [Bond(BondType.TRIPLE)]),
        ("C-C=C", [Bond(BondType.SINGLE), Bond(BondType.DOUBLE)]),
    ],
)
def test_bond(parser, smiles, expected):
    (mol,) = parser.parse(smiles)
    mol.print_table()
    assert bond_values(mol) == expected


@pytest.mark.parametrize(
    "smiles,expected_bonds",
    [
        ("C1C1", 1),
        ("C1CC1", 3),
        # Cubane
        ("C12C3C4C1C5C4C3C25", 12),
        ("CCCB1CCP1CCC", 10),
        ("CCCB=1CCP1CCC", 10),
        ("CCCB1CCP=1CCC", 10),
    ],
)
def test_closure(parser, smiles, expected_bonds):
    (mol,) = parser.parse(smiles)
    mol.print_table()

    assert len(bond_values(mol)) == expected_bonds


@pytest.mark.parametrize(
    "smiles, expected_bonds",
    [
        ("C(C)C", 2),
        ("B=(P)(N)C", 3),
        ("C(C)(C)(C)C", 4),
    ],
)
def test_submolecules(parser, smiles, expected_bonds):
    (mol,) = parser.parse(smiles)
    mol.print_table()

    assert len(bond_values(mol)) == expected_bonds


def test_submolecules_bonds(parser):
    (mol,) = parser.parse("CCCB(CN)(CP)CC")
    mol.print_table()

    assert len(bond_values(mol))


def test_submolecule_cyclo(parser):
    (mol,) = parser.parse("C(=C)(C)C=1CCC1")
    mol.print_table()

    assert len(bond_values(mol)) == 7


def test_submolecule_nested(parser):
    (mol,) = parser.parse("PCCCB(CCCBr(CC)CCN)CCCC")
    mol.print_table()

    assert len(bond_values(mol))


def test_submolecule_cyclo_intersect(parser):
    (mol,) = parser.parse("CCCC(CCB1(CC))CN1CC")
    mol.print_table()
    assert len(bond_values(mol)) == 13


@pytest.mark.parametrize(
    "atom,expected",
    [
        ("[H]", {"symbol": "H"}),
        ("C", {"symbol": "C"}),
        ("[Fe]", {"symbol": "Fe"}),
        ("[22Fe]", {"symbol": "Fe", "isotope": 22}),
        ("[Fe@@]", {"symbol": "Fe", "chirality": "@@"}),
        ("[FeH]", {"symbol": "Fe", "hydrogen": 1}),
        ("[FeH2]", {"symbol": "Fe", "hydrogen": 2}),
        ("[Fe+2]", {"symbol": "Fe", "charge": 2}),
        ("[Fe++]", {"symbol": "Fe", "charge": 2}),
        ("[22Fe+2]", {"symbol": "Fe", "charge": 2, "isotope": 22}),
        (
            "[22Fe@H3+4]",
            {
                "symbol": "Fe",
                "charge": 4,
                "isotope": 22,
                "hydrogen": 3,
                "chirality": "@",
            },
        ),
        (
            "[22Fe@@H3----]",
            {
                "symbol": "Fe",
                "charge": -4,
                "isotope": 22,
                "hydrogen": 3,
                "chirality": "@@",
            },
        ),
        (
            "C",
            {
                "symbol": "C",
                "charge": 0,
                "hydrogen": 4,
            },
        ),
        (
            "[H]",
            {
                "symbol": "H",
                "charge": 0,
                "hydrogen": 1,
            },
        ),
        (
            "N",
            {
                "symbol": "N",
                "charge": 0,
                "hydrogen": 3,
            },
        ),
    ],
)
def test_atom(parser, atom, expected):
    (mol,) = parser.parse(atom)
    mol.print_table()
    atom = mol._atoms[0]

    for key, value in expected.items():
        assert getattr(atom, key) == value


@pytest.mark.parametrize(
    "smiles",
    [
        "C#C#C",
        "C(C)(C)(C)(C)(C)",
        "CCC(#C)CCC",
        "CCC(=N=C)CCC",
        "CC(#N)CC",
    ],
)
def test_assert_valence(parser, smiles):
    with pytest.raises(ValueError):
        (mol,) = parser.parse(smiles)
        mol.print_table()
