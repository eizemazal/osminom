from dataclasses import dataclass
from typing import NotRequired, TypedDict, Unpack

from buftinom.lookup import FunctionalGroup
from buftinom.smileg import Bond, BondSymbol, Molecule
from osminom.atom import Atom


class AtomParams(TypedDict):
    """Params of how to treat the atom during the match process"""

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
    tag: FunctionalGroup


class Matcher:
    """
    This class creates a molecule pattern and,
    given the molecule atom, checks if this pattern is matches at this position.

    """

    def __init__(self, mol: Molecule, atom: AtomParams, tag: FunctionalGroup):
        self.mol = mol
        self.atom = atom
        self.next: list[Matcher] = []
        self.tag = tag

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
        """
        The mach function

        Probe if we can find matcher-defined structure starting from param atom
        """
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

        return GroupMatch(root, atoms, self.tag)


class MatcherBuilder:
    """
    Utility to start building match patterns.
    Attach a tag to the pattern to specify which functional group it matches
    """

    def __init__(self, mol: Molecule, tag: FunctionalGroup):
        self.mol = mol
        self.tag = tag

    def chain(self, *matchers: Matcher):
        m1, *ms = matchers
        root = m1

        for matcher in ms:
            m1 = m1.add(matcher)

        return root

    def atom(self, **atom: Unpack[AtomParams]) -> Matcher:
        return Matcher(self.mol, atom, tag=self.tag)


## Pattern definition
#


def alco_matcher(mol: Molecule):
    match = MatcherBuilder(mol, FunctionalGroup.ALCOHOL)

    matcher = match.chain(
        match.atom(symbol="O"),
        match.atom(by="-", symbol="C", is_root=True),
    )

    return matcher


def acid_matcher(mol: Molecule):
    """Linear match O=C-O"""
    match = MatcherBuilder(mol, FunctionalGroup.CARBOXYLIC_ACID)

    return match.chain(
        match.atom(symbol="O"),
        match.atom(by="=", symbol="C", is_root=True),
        match.atom(by="-", symbol="O"),
    )


def amine_matcher(mol: Molecule):
    match = MatcherBuilder(mol, FunctionalGroup.AMINE)

    return match.chain(
        match.atom(symbol="N"),
        match.atom(by="-", symbol="C", is_root=True),
    )


def acid_matcher2(mol: Molecule):
    """Example of a tree-matcher.
        C = O
       /
      O

    But I recommend to build matchers starting from non-C atoms if possible.
    A tad better performance.
    """
    match = MatcherBuilder(mol, "carboxylic-acid")

    return match.chain(
        match.atom(symbol="C", is_root=True).branch(
            match.atom(by="=", symbol="O"),
            match.atom(by="-", symbol="O"),
        ),
    )


def get_matchers(molecule: Molecule):
    """Return matchers in priority order.

    For examle since acid matcher presceeds alco matcher
    we guarantee that COOH is matching before CO, and selected first (if matched)
    which would cause collision otherwise (because CO always matches if COOH matches)
    """
    return [
        acid_matcher(molecule),
        alco_matcher(molecule),
        amine_matcher(molecule),
    ]
