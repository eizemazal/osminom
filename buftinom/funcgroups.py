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
    is_side_root: NotRequired[bool]


def structcmp(a1: Atom, a2: AtomParams):
    return a1.symbol.lower() == a2["symbol"].lower()


def find(
    mol: Molecule,
    start: Atom,
    via: AtomParams,
    ignore: set[Atom],
) -> Atom | None:
    """
    Looks find neighbour of the `start` atom the matches `via` param
    """
    for a in mol.adj[start]:
        if a in ignore:
            continue
        if structcmp(a, via):
            if "by" not in via:
                return a
            elif mol.bonds[(start, a)] == Bond.make(via["by"]):
                return a


@dataclass(frozen=True, unsafe_hash=True, eq=True)
class GroupMatch:
    root: Atom | None
    side_root: Atom | None
    atoms: set[Atom]
    tag: FunctionalGroup

    def __str__(self):
        return f"group({self.atoms}, {self.tag.name})"

    __repr__ = __str__


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

    def matches(self, start: Atom, *, _visited: set[Atom] = None) -> GroupMatch | None:
        """
        The mach function

        Probe if we can find matcher-defined structure starting from param atom
        """
        _visited = (_visited or set()) | {start}

        if not structcmp(start, self.atom):
            return None

        matches: list[GroupMatch] = []
        for nxt in self.next:
            next_atom = find(self.mol, start, nxt.atom, _visited)
            if next_atom is None:
                return None

            submatch = nxt.matches(next_atom, _visited=_visited)
            if submatch is None:
                return None

            matches.append(submatch)

        root = None
        side_root = None
        atoms: set[Atom] = {start}
        if self.atom.get("is_root"):
            root = start

        if self.atom.get("is_side_root"):
            side_root = start

        for mtch in matches:
            atoms = atoms | mtch.atoms
            if mtch.root is not None:
                if root is not None:
                    raise ValueError(
                        f"Multiple match root found ({root}, {mtch.root}). Use only one atom as root"
                    )

                root = mtch.root

            if mtch.side_root is not None:
                if side_root is not None:
                    raise ValueError(
                        f"Multiple sideroots ({side_root}, {mtch.side_root}) Use only one atom as side_root"
                    )
                side_root = mtch.side_root

        return GroupMatch(root=root, side_root=side_root, atoms=atoms, tag=self.tag)


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
    """C-O"""
    match = MatcherBuilder(mol, FunctionalGroup.ALCOHOL)

    matcher = match.chain(
        match.atom(symbol="O"),
        match.atom(by="-", symbol="C", is_root=True),
    )

    return matcher


def acid_matcher(mol: Molecule):
    """Linear match O=C-O aka R-COOH"""
    match = MatcherBuilder(mol, FunctionalGroup.CARBOXYLIC_ACID)

    return match.chain(
        match.atom(symbol="O"),
        match.atom(by="=", symbol="C", is_root=True),
        match.atom(by="-", symbol="O"),
    )


def amine_matcher(mol: Molecule):
    """C - N"""
    match = MatcherBuilder(mol, FunctionalGroup.AMINE)

    return match.chain(
        match.atom(symbol="N"),
        match.atom(by="-", symbol="C", is_root=True),
    )


def chained_amine_matcher(mol: Molecule):
    """- N -"""
    match = MatcherBuilder(mol, FunctionalGroup.AMINE)

    return match.chain(
        match.atom(symbol="C"),
        match.atom(by="-", symbol="N", is_root=True),
        match.atom(by="-", symbol="C"),
    )


def imine_matcher(mol: Molecule):
    """C = N"""
    match = MatcherBuilder(mol, FunctionalGroup.IMINE)

    return match.chain(
        match.atom(symbol="N"),
        match.atom(by="=", symbol="C", is_root=True),
    )


def nitrile_matcher(mol: Molecule):
    """C # N"""
    match = MatcherBuilder(mol, FunctionalGroup.NITRILE)

    return match.chain(
        match.atom(symbol="N"),
        match.atom(by="#", symbol="C", is_root=True),
    )


def acid_amide_matcher(mol: Molecule):
    """O = C - N"""
    match = MatcherBuilder(mol, FunctionalGroup.AMIDE)

    return match.chain(
        match.atom(symbol="O"),
        match.atom(by="=", symbol="C", is_root=True),
        match.atom(by="-", symbol="N"),
    )


def oxy_matcher(mol: Molecule):
    """- O -"""
    match = MatcherBuilder(mol, FunctionalGroup.OXY)

    return match.chain(
        match.atom(symbol="C"),
        match.atom(by="-", symbol="O", is_root=True),
        match.atom(by="-", symbol="C", is_side_root=True),
    )


def ketone_matcher(mol: Molecule):
    """O = C"""
    match = MatcherBuilder(mol, FunctionalGroup.KETONE)

    return match.chain(
        match.atom(symbol="O"),
        match.atom(by="=", symbol="C", is_root=True),
    )


def ester_matcher(mol: Molecule):
    """O=C-O-C aka RCOOR"""
    match = MatcherBuilder(mol, FunctionalGroup.ESTER)

    return match.chain(
        match.atom(symbol="O"),
        match.atom(by="=", symbol="C", is_root=True),
        match.atom(by="-", symbol="O"),
        match.atom(by="-", symbol="C", is_side_root=True),
    )


def get_matchers(molecule: Molecule):
    """Return matchers in priority order.

    For examle since acid matcher presceeds alco matcher
    we guarantee that COOH is matching before CO, and selected first (if matched)
    which would cause collision otherwise (because CO always matches if COOH matches)
    """
    matchers = [
        ester_matcher,
        acid_amide_matcher,
        acid_matcher,
        oxy_matcher,
        alco_matcher,
        chained_amine_matcher,
        ketone_matcher,
        amine_matcher,
        imine_matcher,
        nitrile_matcher,
    ]

    return [m(molecule) for m in matchers]


# Extra


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
