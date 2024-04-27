from buftinom.smileg import Atom, BondType, Molecule

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


def _bond_valence(mol: Molecule, atom: Atom) -> int:
    bond_valence = 0
    for peer in mol.adj[atom]:
        bond = mol.bonds[(atom, peer)]
        bond_valence += BOND.get(bond.type, 0)

    return bond_valence


def _assert_auto(mol: Molecule, atom: Atom):
    actual_valence = _bond_valence(mol, atom)
    expected = VALENCE[atom.symbol]

    if actual_valence > expected:
        raise ValueError(
            f"Atom {atom} have valence {actual_valence},"
            + f" expected to have no more than {expected}"
        )

    # Calculate missing atom characteristics
    atom.hydrogen = expected - actual_valence
    atom.charge = 0


def _assert_manual(mol: Molecule, atom: Atom):
    bond_valence = _bond_valence(mol, atom)

    expected = VALENCE[atom.symbol]
    hydrogen = getattr(atom, "hydrogen")
    if hydrogen is None:
        hydrogen = 0

    charge = getattr(atom, "charge")
    if charge is None:
        charge = 0

    actual_valence = bond_valence + hydrogen - charge

    # actually incorrect, atom could have more valence than normal,
    # it's just the matter of shell configuration.
    # disable if someone wants to override it manually

    # if actual_valence > expected:
    #     raise ValueError(
    #         f"Atom {atom} have valence {actual_valence}"
    #         + f" ({bond_valence} bonds + {hydrogen}H - {charge} charge), "
    #         + f" expected to have no more than {expected} "
    #     )

    # Calculate missing atom characteristics
    if atom.hydrogen is None:
        atom.hydrogen = expected - actual_valence

    if atom.charge is None:
        atom.charge = 0


def fill_valence(mol: Molecule):
    """
    Asserts correctnes of the atom, based on connection's valence

    Calculates actual valence values for 'auto' atoms
    """
    for atom in mol.atoms:
        if atom.symbol not in VALENCE:
            continue

        if atom.auto:
            _assert_auto(mol, atom)
            continue

        if not atom.auto:
            _assert_manual(mol, atom)
