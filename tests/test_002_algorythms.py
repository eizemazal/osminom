import pytest

from buftinom.iupac import Alogrythms, Chain, MolDecomposition, chainkey
from buftinom.smileg import Atom
from buftinom.smiles_parser import SmilesParser


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

    list = algo.table2list()

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


def test_leaf_cycle():
    algo = algorythms("BCCCCC(C1CN)CC1CP")

    b, n, p = algo.leafs
    assert b.symbol == "B"
    assert n.symbol == "N"
    assert p.symbol == "P"

    d = algo.leaf_distances

    assert d[(b, n)] == 8
    assert d[(b, p)] == 9
    assert d[(n, p)] == 5


def test_max_chain():
    algo = algorythms("BCCCl(CN)CCCP")

    chains = algo.max_chains(algo.all_chains)

    assert len(chains) == 2
    chain = chains[0]

    assert symkey(chain) in {("B", "P"), ("P", "B")}


def test_stipchains():
    algo = algorythms("BCCCl(CN)CCCP")
    chain = algo.max_chain(algo.all_chains)

    subchains = algo.stripchains(algo.all_chains, chain)

    assert len(subchains) == 1

    subchain = subchains[0]

    assert symkey(subchain) in {("Cl", "N"), ("N", "Cl")}


def test_stipchains_many():
    algo = algorythms("BCCCCCCO(C(B)CN)I(C(B)CP)CCCCCCP")
    chain = algo.max_chain(algo.all_chains)

    subchains = algo.stripchains(algo.all_chains, chain)

    assert len(subchains) == 4
    assert set(symkey(sc) for sc in subchains) == {
        ("B", "O"),
        ("N", "O"),
        ("B", "I"),
        ("P", "I"),
    }


def test_stipchains_group():
    algo = algorythms("BCCCCCCO(C(B)CN)I(C(B)CP)CCCCCCP")
    chain = algo.max_chain(algo.all_chains)

    subchains = algo.stripchains(algo.all_chains, chain)
    groupped = algo.group_by_ends(subchains)

    assert len(groupped) == 2
    for _, subchains in groupped.items():
        assert len(subchains) == 2

    assert set(a.symbol for a in groupped.keys()) == {"O", "I"}


@pytest.mark.parametrize(
    "smiles,",
    [
        # They have the same structure, different notation
        "BCCCCCCO(C(B)CN)I(C(B)CP)CCCCCCP",
        "BCCCCCCO(C(B)CN)I(CCCCCCP)C(B)CP",
    ],
)
def test_decompose(smiles):
    algo = algorythms(smiles)

    main = algo.decompose()
    main.print()

    o, i = algo.mol.atom("O"), algo.mol.atom("I")

    assert main.chain[0].symbol == "B"
    assert main.chain[-1].symbol == "P"

    assert len(main.connections) == 2

    o_dec = main.connections[o][0]
    assert o_dec.chain[0].symbol == "C"
    assert o_dec.chain[-1].symbol == "N"

    assert len(o_dec.connections) == 1
    oc_dec = list(o_dec.connections.values())[0][0]
    assert oc_dec.chain[0].symbol == "B"

    i_dec = main.connections[i][0]
    assert i_dec.chain[0].symbol == "C"
    assert i_dec.chain[-1].symbol == "P"

    ic_dec = list(i_dec.connections.values())[0][0]
    assert len(i_dec.connections) == 1
    assert ic_dec.chain[0].symbol == "B"


def test_decompose_multiple_connections():
    algo = algorythms("CCB(C)(C)CC")

    algo.mol.print_table()

    main = algo.decompose()
    main.print()


def test_decompose_multiple_tree_connections():
    algo = algorythms("CCCCB(C(N)C(P)C)(C(C)C(C)C)CCCC")

    main = algo.decompose()
    main.print()


def test_decomps_print():
    algo = algorythms("CCCCCCCCCC")

    a, b, c, d, e, f, g, h, i, j = algo.mol.atoms

    dec = MolDecomposition(
        [a, b, c, d],
        {
            b: [
                MolDecomposition([e, f], {}),
                MolDecomposition([g, h, i], {h: [MolDecomposition([j], {})]}),
            ]
        },
    )

    dec.print()
