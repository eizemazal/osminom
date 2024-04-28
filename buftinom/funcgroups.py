from dataclasses import dataclass
from functools import cached_property
from typing import NotRequired, TypedDict, Unpack

from buftinom.lookup import FunctionalGroup
from buftinom.smileg import Atom, Bond, BondSymbol, Molecule, is_not_carbon


class AtomParams(TypedDict):
    """Params of how to treat the atom during the match process"""

    by: NotRequired[BondSymbol]
    symbol: str
    hydrogen: int
    charge: int

    #
    is_root: NotRequired[bool]
    is_side_root: NotRequired[bool]
    is_terminal: NotRequired[bool]


class GroupParams(TypedDict):
    """Params of the group, how it behaves

    is_symmetric - is it allowed to have root and side_root switch places
        i.e. C-O-CC could be treated as CC-O-C
    """

    is_symmetric: NotRequired[bool]
    is_flex_root: NotRequired[bool]


def structcmp(a1: Atom, a2: AtomParams):
    sym_eq = a1.symbol.lower() == a2["symbol"].lower()
    h_eq = a1.hydrogen == a2.get("hydrogen", a1.hydrogen)
    charge_eq = a1.charge == a2.get("charge", a1.charge)

    return sym_eq and h_eq and charge_eq


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
        if not structcmp(a, via):
            continue

        if "by" not in via:
            return a
        elif mol.bonds[(start, a)] == Bond.make(via["by"]):
            return a


@dataclass(eq=True)
class GroupMatch:
    root: Atom | None
    side_root: Atom | None
    atoms: set[Atom]
    tag: FunctionalGroup
    #
    is_symmetric: bool
    is_flex_root: bool  # can root atom be flexible and connect to the cycle
    #
    root_flexed: Atom | None = None

    @cached_property
    def func_atoms(self) -> set[Atom]:
        return set(filter(is_not_carbon, self.atoms))

    def change_root(self, new_root: Atom):
        if self.root_flexed is not None:
            raise ValueError("Root can only be changed once")
        self.root_flexed = self.root
        self.root = new_root

    def __str__(self):
        return f"group({self.atoms}, {self.tag.name})"

    __repr__ = __str__

    def __hash__(self):
        return hash((self.root, self.side_root, tuple(self.atoms), self.tag))


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
        self.group_params = GroupParams()

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

    def set_params(self, params: GroupParams):
        self.group_params = params

    def matches(self, start: Atom, *, _visited: set[Atom] = None) -> GroupMatch | None:
        """
        The mach function

        Probe if we can find matcher-defined structure starting from param atom
        """
        _visited = (_visited or set()) | {start}

        if not structcmp(start, self.atom):
            return None

        if self.atom.get("is_terminal") and len(self.mol.adj[start]) > 1:
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
                        f"Multiple match root found ({root}, {mtch.root})."
                        + " Use only one atom as root"
                    )

                root = mtch.root

            if mtch.side_root is not None:
                if side_root is not None:
                    raise ValueError(
                        f"Multiple sideroots ({side_root}, {mtch.side_root}) "
                        + "Use only one atom as side_root"
                    )
                side_root = mtch.side_root

        return GroupMatch(
            root=root,
            side_root=side_root,
            atoms=atoms,
            tag=self.tag,
            is_symmetric=self.group_params.get("is_symmetric", False),
            is_flex_root=self.group_params.get("is_flex_root", False),
        )


class MatcherBuilder:
    """
    Utility to start building match patterns.
    Attach a tag to the pattern to specify which functional group it matches
    """

    def __init__(self, mol: Molecule, tag: FunctionalGroup):
        self.mol = mol
        self.tag = tag

    def chain(self, *matchers: Matcher, **group_params: Unpack[GroupParams]) -> Matcher:
        m1, *ms = matchers
        root = m1

        for matcher in ms:
            m1 = m1.add(matcher)

        root.set_params(group_params)
        return root

    def atom(self, **atom: Unpack[AtomParams]) -> Matcher:
        return Matcher(self.mol, atom, tag=self.tag)


## Pattern definition
#


def alco_matcher(mol: Molecule):
    """C-O"""
    match = MatcherBuilder(mol, FunctionalGroup.ALCOHOL)

    matcher = match.chain(
        match.atom(symbol="O", is_terminal=True),
        match.atom(by="-", symbol="C", is_root=True),
    )

    return matcher


def acid_matcher(mol: Molecule):
    """Linear match O=C-O aka R-COOH"""
    match = MatcherBuilder(mol, FunctionalGroup.CARBOXYLIC_ACID)

    return match.chain(
        match.atom(symbol="O", is_terminal=True),
        match.atom(by="=", symbol="C", is_root=True),
        match.atom(by="-", symbol="O", is_terminal=True),
        is_flex_root=True,
    )


def amine_matcher(mol: Molecule):
    """C - N"""
    match = MatcherBuilder(mol, FunctionalGroup.AMINE)

    return match.chain(
        match.atom(symbol="N", is_terminal=True),
        match.atom(by="-", symbol="C", is_root=True),
    )


def imine_matcher(mol: Molecule):
    """C = N"""
    match = MatcherBuilder(mol, FunctionalGroup.IMINE)

    return match.chain(
        match.atom(symbol="N", is_terminal=True),
        match.atom(by="=", symbol="C", is_root=True),
        is_flex_root=True,
    )


def nitrile_matcher(mol: Molecule):
    """C # N"""
    match = MatcherBuilder(mol, FunctionalGroup.NITRILE)

    return match.chain(
        match.atom(symbol="N", is_terminal=True),
        match.atom(by="#", symbol="C", is_root=True),
        is_flex_root=True,
    )


def acid_amide_matcher(mol: Molecule):
    """O = C - N"""
    match = MatcherBuilder(mol, FunctionalGroup.AMIDE)

    return match.chain(
        match.atom(symbol="O", is_terminal=True),
        match.atom(by="=", symbol="C", is_root=True),
        match.atom(by="-", symbol="N", is_terminal=True),
        is_flex_root=True,
    )


def oxy_matcher(mol: Molecule):
    """- O -"""
    match = MatcherBuilder(mol, FunctionalGroup.OXY)

    return match.chain(
        match.atom(symbol="C", is_root=True),
        match.atom(by="-", symbol="O"),
        match.atom(by="-", symbol="C", is_side_root=True),
        is_symmetric=True,
    )


def amino_matcher(mol: Molecule):
    """- N -"""
    match = MatcherBuilder(mol, FunctionalGroup.AMINO)

    return match.chain(
        match.atom(symbol="C", is_root=True),
        match.atom(by="-", symbol="N"),
        match.atom(by="-", symbol="C", is_side_root=True),
        is_symmetric=True,
    )


def ketone_matcher(mol: Molecule):
    """O = C(C)(C)"""
    match = MatcherBuilder(mol, FunctionalGroup.KETONE)

    return match.chain(
        match.atom(symbol="O", is_terminal=True),
        match.atom(by="=", symbol="C", hydrogen=0, is_root=True),
    )


def aldehyde_matcher(mol: Molecule):
    """O = C"""
    match = MatcherBuilder(mol, FunctionalGroup.ALDEHYDE)

    return match.chain(
        match.atom(symbol="O", is_terminal=True),
        match.atom(by="=", symbol="C", hydrogen=1, is_root=True),
        is_flex_root=True,
    )


def brom_matcher(mol: Molecule):
    """Br - C"""
    match = MatcherBuilder(mol, FunctionalGroup.BROM)

    return match.chain(
        match.atom(symbol="Br", is_terminal=True),
        match.atom(by="-", symbol="C", is_root=True),
    )


def chlor_matcher(mol: Molecule):
    """Cl - C"""
    match = MatcherBuilder(mol, FunctionalGroup.CHLOR)

    return match.chain(
        match.atom(symbol="Cl", is_terminal=True),
        match.atom(by="-", symbol="C", is_root=True),
    )


def nitro_matcher(mol: Molecule):
    """R-N(=O)O"""
    match = MatcherBuilder(mol, FunctionalGroup.NITRO)

    return match.chain(
        match.atom(symbol="O", is_terminal=True),
        match.atom(by="=", symbol="N", charge=+1).branch(
            match.atom(by="-", symbol="O", is_terminal=True, charge=-1),
            match.atom(by="-", symbol="C", is_root=True),
        ),
    )


def ester_matcher(mol: Molecule):
    """O=C-O-C aka RCOOR"""
    match = MatcherBuilder(mol, FunctionalGroup.ESTER)

    return match.chain(
        match.atom(symbol="O", is_terminal=True),
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
        nitro_matcher,
        oxy_matcher,
        amino_matcher,
        alco_matcher,
        aldehyde_matcher,
        ketone_matcher,
        amine_matcher,
        imine_matcher,
        nitrile_matcher,
        brom_matcher,
        chlor_matcher,
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
