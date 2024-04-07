from dataclasses import dataclass
from typing import NotRequired, TypedDict, Unpack

from buftinom.smileg import Bond, BondSymbol, Molecule
from osminom.atom import Atom


class AtomParams(TypedDict):
    by: NotRequired[BondSymbol]
    symbol: str
    is_root: NotRequired[bool]


def structcmp(a1: Atom, a2: AtomParams):
    return a1.symbol == a2["symbol"]


def find(mol: Molecule, start: Atom, via: AtomParams) -> Atom | None:
    for a in mol.adj[start]:
        if structcmp(a, via):
            if "by" not in via:
                return a
            elif mol.bonds[(start, a)] == Bond.make(via["by"]):
                return a


@dataclass(frozen=True)
class GroupMatch:
    root: Atom | None
    atoms: set[Atom]


class Matcher:
    def __init__(self, mol: Molecule, atom: AtomParams):
        self.mol = mol
        self.atom = atom
        self.next: list[Matcher] = []

    def then(self, **atom: Unpack[AtomParams]):
        nextm = Matcher(self.mol, atom)
        self.next.append(nextm)
        return nextm

    def add(self, matcher: "Matcher"):
        self.next.append(matcher)
        return matcher

    def branch(self, *matchers: "Matcher"):
        self.next.extend(matchers)
        return self

    def matches(self, start: Atom) -> GroupMatch | None:
        """Test if we can find matcher-defined structure starting from param atom"""
        if not structcmp(start, self.atom):
            return None

        matches: list[GroupMatch] = []
        for nxt in self.next:
            next_atom = find(self.mol, start, nxt.atom)
            if next_atom is None:
                return None

            submatch = nxt.matches(next_atom)
            if submatch is None:
                return None

            matches.append(submatch)

        root = None
        atoms: set[Atom] = {start}
        if self.atom.get("is_root"):
            root = start

        for mtch in matches:
            atoms = atoms | mtch.atoms
            if mtch.root is None:
                continue

            if root is not None:
                raise ValueError(
                    f"Multiple match root found ({root}, {mtch.root}). Use only one atom as root"
                )

            root = mtch.root

        return GroupMatch(root, atoms)


class MatcherBuilder:
    def __init__(self, mol: Molecule):
        self.mol = mol

    def chain(self, *matchers: Matcher):
        m1, *ms = matchers
        root = m1

        for matcher in ms:
            m1 = m1.add(matcher)

        return root

    def atom(self, **atom: Unpack[AtomParams]) -> Matcher:
        return Matcher(self.mol, atom)


def alco_matcher(mol: Molecule):
    match = MatcherBuilder(mol)

    matcher = match.chain(
        match.atom(symbol="O"),
        match.atom(by="-", symbol="C", is_root=True),
    )

    return matcher


def acid_matcher(mol: Molecule):
    match = MatcherBuilder(mol)

    return match.chain(
        match.atom(symbol="O"),
        match.atom(by="=", symbol="C", is_root=True),
        match.atom(by="-", symbol="O"),
    )


def acid_matcher2(mol: Molecule):
    match = MatcherBuilder(mol)

    return match.chain(
        match.atom(symbol="C", is_root=True).branch(
            match.atom(by="=", symbol="O"),
            match.atom(by="-", symbol="O"),
        ),
    )


# def tree_matcher(mol: Molecule, point: Atom):
#     builder = MatcherBuilder(mol)

#     matcher = builder.chain(
#         builder.atom("C"),
#         builder.split(
#             builder.chain(
#                 builder.atom("=", "O"),
#             ),
#             builder.chain(
#                 builder.atom("-", "O"),
#             )
#         )
#     )
