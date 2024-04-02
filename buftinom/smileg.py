from dataclasses import dataclass
from enum import Enum
from typing import Literal


@dataclass(slots=True)
class Atom:
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
        ):
            isotope = f"{self.isotope}" if self.isotope else ""
            chirality = f"@{self.chirality}" if self.chirality else ""
            hydrogens = f"H{self.hydrogen}" if self.hydrogen else ""
            charge = f"{self.charge:+}" if self.charge else ""
            return f"[{isotope}{self.symbol}{chirality}{hydrogens}{charge}]"

        return self.symbol

    __repr__ = __str__

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other


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


@dataclass(slots=True, eq=True, frozen=True, unsafe_hash=True)
class Bond:
    type: BondType

    @classmethod
    def make(cls, s: Literal["-", "=", "#", ":", "/", "\\"]):
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
            raise ValueError(f"Bond `{self._staged_bond}` already staged")

        self._staged_bond = bond
        return self

    def merge(self, other: "Molecule"):
        other.print_table()
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

    @property
    def bonds(self):
        return self._bonds

    def unuque_bonds(self):
        unique_used = set()
        for (a1, a2), b in self._bonds.items():
            if (a1, a2) in unique_used:
                continue
            unique_used.add((a1, a2))
            unique_used.add((a2, a1))
            yield a1, b, a2

    def __str__(self):
        return f"Molecule({len(self._atoms)} atoms, {len(self._bonds)} bonds)"

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

    __repr__ = __str__


def _scenter(obj, w):
    return f"{str(obj):^{w}}"
