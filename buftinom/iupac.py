import heapq
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from functools import cached_property, lru_cache
from typing import Generator, NamedTuple, TypeAlias

from buftinom.lookup import (
    BOND_PRIORITY,
    MULTI_BY_PREFIX,
    MULTI_MULTI_BY_PREFIX,
    ROOT_BY_LENGTH,
    Prefix,
    PrimarySuffix,
    bond_prio_cmp,
)
from buftinom.smileg import Atom, Bond, BondType, Molecule
from buftinom.translate import WordForm

ChainKey: TypeAlias = tuple[Atom, Atom]
Chain: TypeAlias = list[Atom]


def chainkey(chain: Chain) -> ChainKey:
    assert len(chain) > 0, "Empty chain has no key"
    return chain[0], chain[-1]


def reversekey(key: ChainKey) -> ChainKey:
    return key[1], key[0]


def reversechain(chain: Chain) -> Chain:
    return list(reversed(chain))


class FpodReliability(Enum):
    UNKNOWN = "Unknown"
    UNSURE = "Unsure"
    RELIABLE = "100%"


@dataclass
class MolDecomposition:
    chain: Chain
    connections: dict[Atom, list["MolDecomposition"]]
    fpod: FpodReliability = FpodReliability.UNKNOWN

    def print(self, level=0):
        print(self.chain)
        for connection, decomps in self.connections.items():
            for decomp in decomps:
                print("    " * level, connection, ":", end=" ")
                decomp.print(level + 2)

    def subdecompositions(
        self, connector: Atom = None
    ) -> Generator[tuple[Atom | None, "MolDecomposition"], None, None]:
        yield connector, self

        for connector, decomps in self.connections.items():
            for decomp in decomps:
                yield from decomp.subdecompositions(connector)

    def __hash__(self) -> int:
        return id(self)

    def __eq__(self, o: object) -> bool:
        return id(self) == id(o)


@dataclass
class ChainGroup:
    matcher: set[Atom]
    chains: list[Chain]


class Alogrythms:
    def __init__(self, mol: Molecule):
        self.mol = mol
        self.chains: list[Chain] = []

    def table2list(self):
        adj_list = defaultdict(list)

        for a1, a2 in self.mol.bonds.keys():
            adj_list[a1].append(a2)

        return adj_list

    @cached_property
    def adj(self):
        return self.table2list()

    @cached_property
    def leafs(self):
        if len(self.mol.atoms) == 1:
            return [self.mol.atoms[0]]

        return [a for a, lst in self.adj.items() if len(lst) == 1]

    @lru_cache(maxsize=256)
    def distances_from(self, start: Atom):
        """Dijkstra's algorithm"""
        adj = self.adj
        queue = [(0, start)]

        distances = {a: float("inf") for a in self.mol.atoms}
        distances[start] = 0

        predecessors: dict[Atom:Atom] = {a: None for a in self.mol.atoms}

        while queue:
            dist, node = heapq.heappop(queue)

            for neighbor in adj[node]:
                # + 1 here to adjust, for serach with priorities
                new_dist = dist + 1

                if new_dist < distances[neighbor]:
                    distances[neighbor] = new_dist
                    predecessors[neighbor] = node
                    heapq.heappush(queue, (new_dist, neighbor))

        return distances, predecessors

    def unfold_chain(self, predecessors, start, end) -> Chain:
        path = []
        node = end
        while node is not None:
            path.append(node)
            node = predecessors[node]

        path.reverse()
        return path if path[0] == start else []

    @cached_property
    def leaf_distances(self):
        """Calculate distances between all leaf atoms"""
        leafs = self.leafs
        leaf_distances: dict[ChainKey, float] = {}

        for a1 in leafs:
            a1_distances, a1_predecessors = self.distances_from(a1)

            for a2 in leafs:
                if a1 == a2:
                    continue
                self.chains.append(self.unfold_chain(a1_predecessors, a1, a2))
                leaf_distances[(a1, a2)] = a1_distances[a2]

        return leaf_distances

    @cached_property
    def all_chains(self) -> list[Chain]:
        self.leaf_distances
        if len(self.mol.atoms) == 1:
            a = self.mol.atoms[0]
            return [[a]]

        return self.chains

    def max_chains(self, chains: list[Chain]) -> list[Chain]:
        maxlen = len(max(chains, key=len))
        return [chain for chain in chains if len(chain) == maxlen]

    def chain_priority(self, chain: Chain):
        if len(chain) < 2:
            return 0

        prio = 0
        for a1, a2 in zip(chain, chain[1:]):
            bond = self.mol.bonds[(a1, a2)]
            prio += BOND_PRIORITY.get(bond.type, 0)
        return prio

    def max_chain(self, chains: list[Chain]) -> Chain:
        max_chains = self.max_chains(chains)
        #
        max_prio = max(max_chains, key=self.chain_priority)
        #
        return max_prio

    def orient_by_leafs(self, chains: list[Chain]) -> list[Chain]:
        leafs = self.leafs
        oriented: list[Chain] = []

        for chain in chains:
            end = chain[-1]
            if end in leafs:
                chain = reversechain(chain)

            oriented.append(chain)

        return oriented

    def stripchains(self, chains: list[Chain], chain: Chain):
        """Generate all subchains
        All subchains will start with leaf, and end on the input chain
        """
        chain_atoms = set(chain)
        chains = self.orient_by_leafs(chains)
        splits: list[Chain] = []

        existing: set[ChainKey] = {chainkey(chain)}

        for splitting in chains:
            buffer: Chain = []

            for atom in splitting:
                if atom not in chain_atoms:
                    buffer.append(atom)
                    continue

                if not buffer:
                    continue

                buffer.append(atom)
                key = chainkey(buffer)
                if key not in existing:
                    splits.append(list(buffer))
                    existing.add(key)
                    existing.add(reversekey(key))
                    break

                buffer = []

        return splits

    def group_by_ends(self, chains: list[Chain]):
        grouped: dict[Atom, list[Chain]] = defaultdict(list)

        for chain in chains:
            assert len(chain) > 1, chain
            *chain, end = chain

            # reverse to orient the chain from the groupping atom forward
            # no real effect, just easier to understand for humans
            grouped[end].append(reversechain(chain))

        return grouped

    def interconneced(self, chains: list[Chain]) -> list[list[Chain]]:
        if len(chains) == 1:
            return [chains]
        first_chain, *chains = chains
        groups = [ChainGroup(matcher=set(first_chain), chains=[first_chain])]

        for chain in chains:
            chain_atoms = set(chain)
            in_existed_group = False
            for group in groups:
                if chain_atoms & group.matcher:
                    group.chains.append(chain)
                    group.matcher |= set(chain)
                    in_existed_group = True

            if in_existed_group:
                continue

            groups.append(ChainGroup(matcher=set(chain), chains=[chain]))

        return [g.chains for g in groups]

    def decompose(self):
        chains = self.all_chains

        def _decompose(chains: list[Chain]):
            connections: dict[Atom, list[Chain]] = defaultdict(list)

            max_chain = self.max_chain(chains)
            subchains = self.stripchains(chains, max_chain)
            groupped = self.group_by_ends(subchains)

            for connection, groupped_chains in groupped.items():
                for chaingroup in self.interconneced(groupped_chains):
                    connections[connection].append(_decompose(chaingroup))

            return MolDecomposition(tuple(max_chain), connections)

        return _decompose(chains)


@dataclass
class AtomFeatures:
    chain_index: int
    atom: Atom
    bond_ahead: Bond
    subchains: list[MolDecomposition]


class Synt(NamedTuple):
    """SYNtax uniT"""

    index: int
    value: WordForm

    def __str__(self):
        return f"{self.index}-{self.value}"

    __repr__ = __str__


class Subn(NamedTuple):
    """SUB Name"""

    index: int
    value: "IupacName"

    def __str__(self):
        return f"{self.index}-{self.value}"

    __repr__ = __str__


def s(obj):
    if not obj:
        return ""
    return str(obj)


@dataclass(frozen=True, slots=True, kw_only=True)
class IupacName:
    prefixes: list[Subn] = field(default_factory=list)
    infix: WordForm = None
    root_word: WordForm
    prime_suffixes: list[Synt]
    sub_suffix: WordForm = None
    func_suffix: WordForm = None

    ref: MolDecomposition

    def __repr__(self):
        return "".join(
            map(
                s,
                [
                    self.prefixes,
                    self.infix,
                    self.root_word,
                    self.prime_suffixes,
                    self.sub_suffix,
                    self.func_suffix,
                ],
            )
        )


def suffixes2str(suffixes: list[Synt], sub_suffix: WordForm):
    res: list[str] = []
    ssuff: list[Synt] = sorted(suffixes, key=lambda s: s.value)
    for idx, name in ssuff:
        if idx > 1:
            res.extend(["-", str(idx), "-", name.short])
        else:
            res.append(name.short)

    if sub_suffix:
        res.extend(sub_suffix.norm)
    else:
        res[-1] = ssuff[-1].value.norm

    return "".join(res)


def revsubnames(chain_len: int, subnames: list[tuple[int, str]]):
    for idx, name in subnames:
        yield chain_len - idx + 1, name


def fpod2(chain_len: int, subnames: list[tuple[int, str]]):
    """First point of difference second round to find weher we need to reverse chain
    with a new information about the names of the subchains.

    By this time we know that we really need this info
    to resolve possible collisions (see Iupac.fpod)
    """
    idmap = defaultdict(list)
    for idx, name in subnames:
        idmap[idx].append(name)

    for idx in idmap.keys():
        idmap[idx] = tuple(sorted(idmap[idx]))

    for idx in sorted(idmap.keys()):
        rev_idx = chain_len - idx + 1

        cur = idmap[idx]
        rev = idmap.get(rev_idx, None)
        if not rev:
            return subnames

        if cur < rev:
            # all good, alphabetical order
            return subnames

        if cur > rev:
            # 🔥 Collision, reverse indexing have better alphabetical order
            return list(revsubnames(chain_len, subnames))

    return subnames


def preffixes2str(iupac: IupacName):
    preffixes = iupac.prefixes
    if not preffixes:
        return None

    #
    subnames = [(idx, iupac2str(subn)) for idx, subn in preffixes]
    #
    if iupac.ref.fpod != FpodReliability.RELIABLE:
        subnames = fpod2(len(iupac.ref.chain), subnames)

    groups: dict[str, list[int]] = defaultdict(list)

    for idx, name in subnames:
        groups[name].append(idx)

    ordered = sorted(list(groups.items()), key=lambda g: g[0])

    res: list[str] = []

    multiselector = MULTI_BY_PREFIX

    for name, indexes in ordered:
        subname: list[str] = []

        name_complex = name[0].isdigit()
        if name_complex:
            multiselector = MULTI_MULTI_BY_PREFIX

        if indexes:
            indexes = list(map(str, sorted(indexes)))
            subname.extend(
                [
                    ",".join(indexes),
                    "-",
                    multiselector.get(len(indexes)).value.norm,
                ]
            )

        if name_complex:
            subname.append("(")

        subname.append(name)

        if name_complex:
            subname.append(")")

        res.append("".join(subname))

    return "-".join(res)


def iupac2str(iupac: IupacName):
    preffix = preffixes2str(iupac)
    infix = iupac.infix
    root = iupac.root_word.norm
    suffix = suffixes2str(iupac.prime_suffixes, iupac.sub_suffix)
    fsuffix = iupac.func_suffix

    return "".join(map(s, [preffix, infix, root, suffix, fsuffix]))


class Iupac:
    def __init__(self, mol: Molecule):
        self.mol = mol

    @cached_property
    def decomposition(self):
        return Alogrythms(self.mol).decompose()

    def features(self, decomposition: MolDecomposition):
        """Excluding last atom"""
        chain = decomposition.chain

        for i, (a1, a2) in enumerate(zip(chain, chain[1:]), start=1):
            bond = self.mol.bonds[(a1, a2)]

            yield AtomFeatures(i, a1, bond, decomposition.connections.get(a1, []))

    def fpod(self, decomposition: MolDecomposition):
        straight = decomposition
        reverse = MolDecomposition(
            tuple(reversed(straight.chain)), straight.connections
        )

        straight_features = list(self.features(straight))
        reverse_features = list(self.features(reverse))

        nsubchains = 0
        nbonds = 0
        rel = FpodReliability.RELIABLE

        for fstraight, freversed in zip(straight_features, reverse_features):
            sub_strait = len(fstraight.subchains)
            sub_revers = len(freversed.subchains)

            # First subchain
            subchain_cmp = sub_strait - sub_revers
            if subchain_cmp > 0:
                return straight, straight_features, rel
            if subchain_cmp < 0:
                return reverse, reverse_features, rel

            # First posin of unsaturation
            bond_cmp = bond_prio_cmp(fstraight.bond_ahead, freversed.bond_ahead)
            if bond_cmp > 0:
                return straight, straight_features, rel
            if bond_cmp < 0:
                return reverse, reverse_features, rel

            nsubchains += sub_strait
            if fstraight.bond_ahead.type != BondType.SINGLE:
                nbonds += 1

        if nbonds + nsubchains > 0:
            rel = FpodReliability.UNSURE

        return straight, straight_features, rel

    def chain_name(decomp: MolDecomposition, features: list[AtomFeatures]):
        pass

    def suffixes_by_features(self, features: list[AtomFeatures], *, primary: bool):
        result = []

        for feature in features:
            if feature.bond_ahead.type == BondType.DOUBLE:
                result.append(Synt(feature.chain_index, PrimarySuffix.ENE.value))

            if feature.bond_ahead.type == BondType.TRIPLE:
                result.append(Synt(feature.chain_index, PrimarySuffix.YNE.value))

        if not result and primary:
            result.append(Synt(1, PrimarySuffix.ANE.value))

        return result

    def preffixes_by_features(self, features: list[AtomFeatures]):
        result = []

        for feature in features:
            for subchain in feature.subchains:
                #
                subname = self.decompose_name(subchain, primary=False)
                #
                result.append(Subn(feature.chain_index, subname))
        return result

    def decompose_name(self, decomposition: MolDecomposition, *, primary: bool = True):
        decomp, features, reliability = self.fpod(decomposition)
        decomp.fpod = reliability

        root = ROOT_BY_LENGTH[len(decomp.chain)].value
        suffix = self.suffixes_by_features(features, primary=primary)
        subsuff = Prefix.YL.value if not primary else None
        preffix = self.preffixes_by_features(features)

        return IupacName(
            prefixes=preffix,
            infix=None,
            root_word=root,
            prime_suffixes=suffix,
            sub_suffix=subsuff,
            func_suffix=None,
            ref=decomp,
        )
