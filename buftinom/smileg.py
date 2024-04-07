from ast import TypeAlias
from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from typing import Literal

DEBUG_ATOMS = False


def debug_atoms(debug: bool):
    """Wether to print atom id in __str__ and __repr__"""
    global DEBUG_ATOMS
    DEBUG_ATOMS = True


@dataclass(slots=True, frozen=True)
class Atom:
    id: int
    symbol: str
    isotope: int | None = None
    chirality: str | None = None
    hydrogen: int | None = None
    charge: int = 0

    def __str__(self):
        if (
            self.isotope is not None
            or self.hydrogen is not None
            or self.chirality is not None
            or self.charge != 0
        ):
            isotope = f"{self.isotope}" if self.isotope else ""
            chirality = f"@{self.chirality}" if self.chirality else ""
            hydrogens = f"H{self.hydrogen}" if self.hydrogen else ""
            charge = f"{self.charge:+}" if self.charge else ""
            if self.charge and abs(self.charge) == 1:
                charge = charge.replace("1", "")
            return f"[{isotope}{self.symbol}{chirality}{hydrogens}{charge}]"

        if DEBUG_ATOMS:
            return f"{self.symbol}'{self.id}"

        return self.symbol

    __repr__ = __str__

    def __hash__(self):
        return id(self.id)

    def __eq__(self, other):
        if isinstance(other, Atom):
            return self.id == other.id
        return False

    def __lt__(self, other):
        return self.id < other.id


class BondType(Enum):
    SINGLE = "-"
    DOUBLE = "="
    TRIPLE = "#"
    AROMATIC = ":"
    CIS = "/"
    TRANS = "\\"

    def __str__(self):
        return self.value

    __repr__ = __str__


BondSymbol: TypeAlias = Literal["-", "=", "#", ":", "/", "\\"]


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
        self._staged_bond: Bond | None = None
        self._bonds: dict[tuple[Atom, Atom] : Bond] = {}
        self._closures: dict[str : tuple[Bond, Atom]] = {}

    @property
    def bonds(self) -> dict[tuple[Atom, Atom], Bond]:
        return self._bonds

    @property
    def atoms(self) -> list[Atom]:
        return self._atoms

    def atom(self, symbol: str):
        """find first atom by symbol, for test purposes"""
        return [a for a in self._atoms if a.symbol == symbol][0]

    def table2list(self):
        adj_list: dict[Atom, list[Atom]] = {}

        for a in self._atoms:
            adj_list[a] = []

        for a1, a2 in self.bonds.keys():
            adj_list[a1].append(a2)

        return adj_list

    @cached_property
    def adj(self):
        return self.table2list()

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


class MoleculeConstructor(Molecule):
    def __init__(self):
        super().__init__()
        self.connect_bond = Bond.make("-")

    def add_atom(self, atom: Atom):
        assert atom not in self._atoms, f"Atom {atom} already exists"
        if self._atoms:
            raise ValueError("It is not the first atom, use bind_atom instead")

        self._atoms.append(atom)
        self._bind_target = atom
        return self

    def bind_atom(self, atom: Atom):
        bond = self._pop_bond()
        self._atoms.append(atom)
        self._add_bond(self._bind_target, bond, atom)
        self._bind_target = atom

        return self

    def bind_closure(self, closure_id: str):
        bond = self._pop_bond()
        return self._add_closure(closure_id, bond, self._bind_target)

    def _add_closure(self, closure_id: str, bond: Bond, atom: Atom):
        if closure_id not in self._closures:
            self._closures[closure_id] = (bond, atom)
            return self

        first_bond, target = self._closures.pop(closure_id)
        selected_bond = self._select_bond(first_bond, bond)
        self._add_bond(atom, selected_bond, target)
        return self

    def push_bond(self, bond: Bond):
        if self._staged_bond is not None:
            raise ValueError(
                f"Invalid syntax: multiple bonds `{self._staged_bond}, {bond}`"
            )

        self._staged_bond = bond
        return self

    def merge(self, other: "MoleculeConstructor"):
        cbond = self._pop_bond()
        incoming_bond = other.connect_bond
        connectee = other._atoms[0]

        self._atoms.extend(other._atoms)
        self._bonds.update(other._bonds)

        self._add_bond(
            self._bind_target, self._select_bond(cbond, incoming_bond), connectee
        )

        for closure_id, (bond, target) in other._closures.items():
            self._add_closure(closure_id, bond, target)

        return self

    def _pop_bond(self):
        if self._staged_bond is None:
            return Bond.make("-")

        bond = self._staged_bond
        self._staged_bond = None
        return bond

    def _select_bond(self, first_bond: Bond, new_bond: Bond):
        if first_bond.type == BondType.SINGLE:
            return new_bond

        if new_bond.type == BondType.SINGLE:
            return first_bond

        # first_bond have priority in ambiguious choices, unless single
        # Syntax error, actually
        return first_bond

    def _add_bond(self, atom1: Atom, bond: Bond, atom2: Atom):
        if (atom1, atom2) in self._bonds:
            if self._bonds[(atom1, atom2)] == bond:
                return self
            raise ValueError(f"Bond {atom1}{bond}{atom2} already exists")

        self._bonds[(atom1, atom2)] = bond
        self._bonds[(atom2, atom1)] = bond
        return self

    def unuque_bonds(self):
        unique_used = set()
        for (a1, a2), b in self._bonds.items():
            if (a1, a2) in unique_used:
                continue
            unique_used.add((a1, a2))
            unique_used.add((a2, a1))
            yield a1, b, a2


def _scenter(obj, w):
    return f"{str(obj):^{w}}"
