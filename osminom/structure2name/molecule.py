# Graph representation of the molecule
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from typing import Literal
from osminom.common.atom import Atom
from osminom.common.ast import Ast, AstNode
from osminom.common.smiles_parser import SmilesParser


ELEMENT_VALENCE = {
    "C": 4,
    "H": 1,
    "O": 2,
    "N": 3,
}


DEBUG_ATOMS = False


def debug_atoms(debug: bool):
    """Wether to print atom id in __str__ and __repr__"""
    global DEBUG_ATOMS
    DEBUG_ATOMS = True


def s(val, template):
    if not val:
        return ""

    return template % val


type BondSymbol = Literal["-", "=", "#", ":", "/", "\\"]


class BondType(Enum):
    SINGLE = "-"
    DOUBLE = "="
    TRIPLE = "#"
    AROMATIC = ":"
    CIS = "/"
    TRANS = "\\"

    def __str__(self):
        return self.value

    @staticmethod
    def valence(bs: BondSymbol):
        return {
            BondType.SINGLE: 1,
            BondType.DOUBLE: 2,
            BondType.TRIPLE: 3,
        }.get(bs, 1)

    __repr__ = __str__


@dataclass(slots=True, eq=True, frozen=True, unsafe_hash=True)
class Bond:
    type: BondType

    @classmethod
    def make(cls, s: BondSymbol):
        return cls(BondType(s))

    def __str__(self):
        return self.type.value

    __repr__ = __str__


class Molecule:
    def __init__(self):
        self._atoms: list[Atom] = []
        self._bind_target: Atom | None = None
        self._upbond: str | None = None
        self._staged_bond: Bond | None = None
        self._bonds: dict[tuple[Atom, Atom], Bond] = {}
        self._closures: dict[str, tuple[Bond, Atom]] = {}

    @property
    def bonds(self) -> dict[tuple[Atom, Atom], Bond]:
        return self._bonds

    @property
    def atoms(self) -> list[Atom]:
        return self._atoms

    def atom(self, symbol: str):
        """find first atom by symbol, for test purposes"""
        return [a for a in self._atoms if a.symbol == symbol][0]

    def _adjacent_list(self):
        adj_list: dict[Atom, list[Atom]] = {}

        for a in self._atoms:
            adj_list[a] = []

        for a1, a2 in self.bonds.keys():
            adj_list[a1].append(a2)

        return adj_list

    @cached_property
    def adj(self):
        return self._adjacent_list()

    @cached_property
    def adj_set(self):
        return {a: set(conns) for a, conns in self.adj.items()}

    def print_table(self):
        w = max(len(str(a)) for a in self._atoms) + 2

        print(_scenter("\\", w), end="")
        for a in self._atoms:
            print(_scenter(a, w), end="")
        print()

        for a1 in self._atoms:
            print(_scenter(a1, w), end="")
            for a2 in self._atoms:
                bond = self._bonds.get((a1, a2), " ")
                print(_scenter(bond, w), end="")
            print()
        print()

    def __str__(self):
        return f"Molecule({len(self._atoms)} atoms, {len(self._bonds) // 2} bonds)"

    __repr__ = __str__

    def bind(self, atom1: Atom, atom2: Atom, bond_str: str) -> None:
        self._bonds[(atom1, atom2)] = self._bonds[(atom2, atom1)] = Bond.make(bond_str)

    @staticmethod
    def accumulate_fn(this_node: AstNode, downwards: list[Molecule]) -> Molecule:
        mol = Molecule()
        mol._atoms.append(this_node.atom)
        mol._upbond = this_node.upbond
        for downward in downwards:
            mol._atoms.extend(downward._atoms)
            mol._bonds.update(downward._bonds)
            mol.bind(this_node.atom, downward._atoms[0], downward._upbond)
        return mol

    @staticmethod
    def from_ast(ast: Ast) -> Molecule:
        product = ast.root.accumulate(Molecule.accumulate_fn)
        closure_bonds: dict[str, str] = {}
        closure_atoms: dict[str, list[Atom]] = {}
        for node in ast.root.traverse():
            for id, bond_str in node.atom.closures.items():
                existing = closure_bonds.get(id, bond_str)
                closure_bonds[id] = existing if existing != "-" else bond_str
                closure_atoms.setdefault(id, []).append(node.atom)
        if set(closure_atoms.keys()) != set(closure_bonds.keys()):
            raise ValueError("Mismatch in closure bonds and atoms - this cannot be")
        if any([len(cla) != 2 for cla in closure_atoms.values()]):
            raise ValueError("Unpaired atoms in closure list")
        for idx, bond_str in closure_bonds.items():
            product.bind(closure_atoms[idx][0], closure_atoms[idx][1], bond_str)
        product.add_hydrogens()
        return product

    @staticmethod
    def from_smiles(smiles: str) -> Molecule:
        ast = SmilesParser().parse(smiles)
        return Molecule.from_ast(ast)

    def add_hydrogens(self):
        for atom in self._atoms:
            if atom.explicit_protons:
                continue
            atom_valence = ELEMENT_VALENCE.get(atom.symbol)
            if atom_valence is None:
                continue
            if atom.charge == 1 and atom.symbol.upper() in ("N", "P"):
                # onium salts like NH4+
                atom_valence += 1
            elif atom.charge in (-1, -2) and atom.symbol in ("Si"):
                # ate salts like [SiF6]2-
                atom_valence -= atom.charge
            else:
                atom_valence -= abs(atom.charge)

            covalent_bond_count = sum(
                [
                    BondType.valence(bond.type)
                    for atoms, bond in self._bonds.items()
                    if atoms[0] is atom
                ]
            )
            nprotons = atom_valence - covalent_bond_count
            if nprotons < 0:
                raise ValueError(
                    f"Unexpected valence {covalent_bond_count} for element {atom.symbol}; specify explicit protons to avoid this error."
                )
            atom.nprotons = nprotons

    def unique_bonds(self):
        unique_used = set()
        for (a1, a2), b in self._bonds.items():
            if (a1, a2) in unique_used:
                continue
            unique_used.add((a1, a2))
            unique_used.add((a2, a1))
            yield a1, b, a2


def _scenter(obj, w):
    return f"{str(obj):^{w}}"
