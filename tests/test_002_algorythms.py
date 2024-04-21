import pytest

from buftinom.algorythms import Alogrythms, Chain, chainkey
from buftinom.smileg import Atom, debug_atoms
from buftinom.smiles_parser import SmilesParser

debug_atoms(True)


def algorythms(smiles: str):
    return Alogrythms(SmilesParser(debug=True).parse(smiles)[0])


def a2s(*atoms: Atom):
    return tuple([a.symbol for a in atoms])


def symkey(chain: Chain):
    a, b = chainkey(chain)
    return a.symbol, b.symbol


@pytest.mark.parametrize(
    "smiles,expected_length",
    [
        ("C", 1),
        ("CC", 2),
        ("CCCCCCC", 2),
        ("CCC(CC)CCCCC", 3),
        ("CCC(CC)(CC)CCCC", 4),
    ],
)
def test_leafs(smiles, expected_length):
    assert len(algorythms(smiles).leafs) == expected_length


def test_table2list():
    algo = algorythms("CCC")
    a1 = algo.mol._atoms[0]
    a2 = algo.mol._atoms[1]
    a3 = algo.mol._atoms[2]

    list = algo.mol.table2list()

    assert len(list[a1]) == 1
    assert len(list[a2]) == 2
    assert len(list[a3]) == 1


def test_leaf_distances():
    algo = algorythms("CCCCC")
    a1, a5 = algo.leafs

    d = algo.leaf_distances

    assert d[(a1, a5)] == 4


def test_leaf_distances_tree():
    algo = algorythms("CCC(CC)(CC)CC")

    d = algo.leaf_distances

    for a1, a2 in zip(algo.leafs, algo.leafs):
        if a1 == a2:
            continue
        assert d[(a1, a2)] == 4


@pytest.mark.parametrize(
    "smiles, expected_len",
    [
        ("CCCC(CC)CCCC", 8),
        # Not 7 bc O is the functional group
        ("CCC(C)C(C(C)CC)C(O)C", 6),
    ],
)
def test_max_chain(smiles, expected_len):
    algo = algorythms(smiles)

    chain = algo.max_chain(algo.all_chains)

    assert len(chain) == expected_len


def test_stipchains():
    algo = algorythms("CCCC(CC)CCCC")
    chain = algo.max_chain(algo.all_chains)

    subchains = algo.stripchains(algo.all_chains, chain)

    assert len(subchains) == 1


def test_stipchains_many():
    algo = algorythms("CCCCCCCC(C(C)CC)C(C(C)CC)CCCCCCC")
    chain = algo.max_chain(algo.all_chains)

    subchains = algo.stripchains(algo.all_chains, chain)

    assert len(subchains) == 6


@pytest.mark.parametrize(
    "smiles,",
    [
        # They have the same structure, different notation
        "CCCCCCCC(C(C)CC)C(C(C)CC)CCCCCCC",
        "CCCCCCCC(C(C)CC)C(CCCCCCC)C(C)CC",
    ],
)
def test_decompose(smiles):
    algo = algorythms(smiles)

    main = algo.decompose()
    main.print()

    assert len(main.connections) == 2


def test_decompose_multiple_connections():
    algo = algorythms("CCC(C)(C)CC")

    algo.mol.print_table()

    main = algo.decompose()
    main.print()


def test_decompose_multiple_tree_connections():
    algo = algorythms("CCCCC(C(C)C(C)C)(C(C)C(C)C)CCCC")

    main = algo.decompose()
    main.print()


@pytest.mark.parametrize(
    "smiles, start, expected_len",
    [
        ("CCCCCCCCCCCC", "C", 0),
        ("C", "C", 0),
        ("C1C1", "C", 0),
        ("C1CC1", "C", 3),
        ("C1CCCCC1", "C", 6),
        # with subchains
        ("CC1CCCCC1", "C", 6),
        ("CCCC1CCCCC1", "C", 6),
        ("CC1CCCCC1CC", "C", 6),
        ("CC1CC(CCC)CCC1CC", "C", 6),
        # inside the circle
        ("CCCC1CBCCC1", "B", 6),
    ],
)
def test_cycle(smiles, start, expected_len):
    debug_atoms(True)
    algo = algorythms(smiles)
    cycle = algo.shortest_cycle(algo.mol.adj, algo.mol.atom(start))
    print(cycle)
    assert len(cycle) == expected_len, cycle


@pytest.mark.parametrize(
    "smiles, expected_len",
    [
        ("CC1CCCCC1", 1),
        ("C1CCCCC1C1CCCCC1", 2),
        # just plot it :D
        ("C1CCB3P1CCCCCCCCCCCCCCCCCCCCCCCI2CCCCN2CCCCCCCCCCCCCCCCCCCCCCCC3", 3),
    ],
)
def test_cycles(smiles, expected_len):
    debug_atoms(True)
    algo = algorythms(smiles)
    cycles = algo.cycles

    assert len(cycles) == expected_len


@pytest.mark.parametrize(
    "smiles",
    [
        "CC1CCCCC1",
        "CC1C(C)CC(C)CC1",
        # "C1CCCCC1C1CCCCC1",
    ],
)
def test_decompose_cycles(smiles):
    debug_atoms(True)
    algo = algorythms(smiles)
    dec = algo.decompose()

    dec.print()


@pytest.mark.parametrize(
    "smiles, count",
    [
        ("CCC(=O)O", 1),
        ("CCCO", 1),
        ("CC(O)C", 1),
        ("CCCN", 1),
    ],
)
def test_functional_groups(smiles, count):
    debug_atoms(True)
    algo = algorythms(smiles)
    groups = algo.functional_groups

    for g in groups.items():
        print(g)
    assert len(groups) == count


@pytest.mark.parametrize(
    "smiles, count",
    [
        ("CCC(=O)O", 3),
        ("CCCO", 3),
        ("CC(O)C", 3),
        ("CCCN", 3),
    ],
)
def test_one_functional_group_decompostion(smiles, count):
    debug_atoms(True)
    algo = algorythms(smiles)
    dec = algo.decompose()

    dec.print()

    assert len(dec.chain) == count
    assert len(dec.functional_groups) == 1
