from buftinom.smileg import BondType, Molecule

VALENCE = {
    "C": 4,
    "H": 1,
    "O": 2,
    "N": 3,
}


BOND = {
    BondType.SINGLE: 1,
    BondType.DOUBLE: 2,
    BondType.TRIPLE: 3,
}


def assert_valence(mol: Molecule):
    adj = mol.adj

    for atom in mol.atoms:
        if atom.symbol not in VALENCE:
            continue

        actual_valence = 0

        for peer in adj[atom]:
            bond = mol.bonds[(atom, peer)]

            actual_valence += BOND.get(bond.type, 0)

        expected = VALENCE[atom.symbol]

        if actual_valence > expected:
            raise ValueError(
                f"Atom {atom} have valence {actual_valence},"
                + f" expected to have no more than {expected}"
            )
