import heapq
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property, lru_cache
from typing import Generator, NotRequired, TypeAlias, TypedDict, Unpack

from buftinom.lookup import ROOT_BY_LENGTH, Prefix, PrimarySuffix, bond_prio_cmp
from buftinom.smileg import Atom, Bond, BondType, Molecule, MoleculeConstructor

ChainKey: TypeAlias = tuple[Atom, Atom]
Chain: TypeAlias = list[Atom]


def chainkey(chain: Chain) -> ChainKey:
    assert len(chain) > 0, "Empty chain has no key"
    return chain[0], chain[-1]


def reversekey(key: ChainKey) -> ChainKey:
    return key[1], key[0]


def reversechain(chain: Chain) -> Chain:
    return list(reversed(chain))


@dataclass(frozen=True)
class MolDecomposition:
    chain: Chain
    connections: dict[Atom, list["MolDecomposition"]]

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

    def _unfold_chain(self, predecessors, start, end) -> Chain:
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
                self.chains.append(self._unfold_chain(a1_predecessors, a1, a2))
                leaf_distances[(a1, a2)] = a1_distances[a2]

        return leaf_distances

    @cached_property
    def deduped_chains(self) -> list[Chain]:
        self.leaf_distances
        if len(self.mol.atoms) == 1:
            a = self.mol.atoms[0]
            return [[a]]

        deduped: list[Chain] = []

        existing_chains: set[ChainKey] = set()

        for chain in self.chains:
            key = chainkey(chain)
            if key in existing_chains:
                continue

            existing_chains.add(key)
            existing_chains.add(reversekey(key))

            deduped.append(chain)

        return deduped

    def max_chains(self, chains: list[Chain]) -> list[Chain]:
        maxlen = len(max(chains, key=len))
        return [chain for chain in chains if len(chain) == maxlen]

    def max_chain(self, chains: list[Chain]) -> Chain:
        max_chains = self.max_chains(chains)
        # assert len(max_chains) == 1, "Not implemented for multiple max chains"
        # first is ok for not but there are priorities
        max_chain = max_chains[0]
        return max_chain

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

        existing: set[ChainKey] = set()

        for splitting in chains:
            buffer: Chain = []

            for atom in splitting:
                if atom not in chain_atoms:
                    buffer.append(atom)
                    continue

                if not buffer:
                    continue

                key = chainkey(buffer)
                if key not in existing:
                    buffer.append(atom)
                    splits.append(list(buffer))
                    existing.add(key)
                    existing.add(reversekey(key))

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
        chains = self.deduped_chains

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
class ChainGroup:
    matcher: set[Atom]
    chains: list[Chain]


class Features:
    def __init__(self, mol: MoleculeConstructor, chain: list[Atom]):
        self.mol = mol
        self.chain = chain

    def length(self):
        return len(self.chain)

    def fpod(self):
        """first point of difference"""
        for x in zip(
            self.unusual_bonds(reverse=False), self.unusual_bonds(reverse=True)
        ):
            (i, (_, a, _)), (j, (_, b, _)) = x
            if i < j:
                break
            if i > j:
                self.chain = list(reversed(self.chain))
                break
            if i == j:
                bond_prio = bond_prio_cmp(a, b)
                if bond_prio > 0:
                    break
                if bond_prio < 0:
                    self.chain = list(reversed(self.chain))
                    break
                break

    def bonds(self, reverse):
        if len(self.chain) == 1:
            return

        chain = self.chain
        if reverse:
            chain = list(reversed(chain))

        for i, (a1, a2) in enumerate(zip(chain, chain[1:])):
            bond = self.mol.bonds.get((a1, a2))
            if not bond:
                raise ValueError(f"Missing bond between {a1} and {a2}")

            yield i, (a1, bond, a2)

    def unusual_bonds(self, reverse=False):
        for i, (a1, bond, a2) in self.bonds(reverse):
            if bond.type != BondType.SINGLE:
                yield i, (a1, bond, a2)


class NamingOptionsDict(TypedDict):
    primary_chain: NotRequired[bool]


@dataclass(frozen=True)
class NamingOptions:
    primary_chain: bool = False


class ChainNamer:
    def __init__(
        self, mol: Molecule, dec: MolDecomposition, **options: Unpack[NamingOptionsDict]
    ):
        self.mol = mol
        self.dec = dec
        self.options = NamingOptions(**options)

    def suffix_by_bonds(self, bonds: list[tuple[int, tuple[Atom, BondType, Atom]]]):
        if not bonds:
            yield (0, PrimarySuffix.ANE)

        for i, (a1, bond, a2) in bonds:
            match bond.type:
                case BondType.SINGLE:
                    yield (i, PrimarySuffix.ANE)
                case BondType.DOUBLE:
                    yield (i, PrimarySuffix.ENE)
                case BondType.TRIPLE:
                    yield (i, PrimarySuffix.YNE)

    def primary_suffixes(self, primary_suffixes: list[tuple[int, PrimarySuffix]]):
        result = []

        for i, (pos, suffix) in enumerate(primary_suffixes):
            form = suffix.value.short
            if i == len(primary_suffixes) - 1:
                if self.options.primary_chain:
                    form = suffix.value.norm
                else:
                    form = Prefix.YL.value.norm

            if pos == 0:
                result.append(form)
                continue

            result.extend([str(pos + 1), form])

        return result

    def simple_chain_name(self, chain: list[Atom]):
        f = Features(self.mol, chain)
        f.fpod()

        root = ROOT_BY_LENGTH[f.length()]

        unusual_bonds = list(f.unusual_bonds())
        primary_suffixes = list(self.suffix_by_bonds(unusual_bonds))
        primary_suffixes = sorted(primary_suffixes, key=lambda x: x[1].value.norm)

        psuffixes = self.primary_suffixes(primary_suffixes)
        joined = "-".join(psuffixes)
        link = "-" if psuffixes[0].isnumeric() else ""

        return f"{root.value.norm}{link}{joined}"


@dataclass
class AtomFeatures:
    chain_index: int
    atom: Atom
    bond_ahead: Bond
    subchain: MolDecomposition | None


class Iupac:
    def __init__(self, mol: Molecule):
        self.mol = mol

    @cached_property
    def decomposition(self):
        return Alogrythms(self.mol).decompose()

    def subchain_simple_names(self):
        decomposition = self.decomposition

        # decomp : (ConnectorAtom, str)
        # ConnectorAtom == None for primary chain

        for connector, dec in decomposition.subdecompositions():
            namer = ChainNamer(
                self.mol,
                dec,
                primary_chain=dec == decomposition,
            )

            yield dec, connector, namer.simple_chain_name(dec.chain)

    def features(self, decomposition: MolDecomposition):
        chain = decomposition.chain

        for i, (a1, a2) in enumerate(zip(chain, chain[1:]), start=1):
            bond = self.mol.bonds[(a1, a2)]

            yield AtomFeatures(i, a1, bond, decomposition.connections.get(a2))

    def fpod(self, decomposition: MolDecomposition):
        straight = decomposition
        reverse = MolDecomposition(
            tuple(reversed(straight.chain)), straight.connections
        )

        straight_features = list(self.features(straight))
        reverse_features = list(self.features(reverse))

        for fstraight, freversed in zip(straight_features, reverse_features):
            substrait = fstraight.subchain is not None
            subrevers = freversed.subchain is not None

            bond_cmp = bond_prio_cmp(fstraight.bond_ahead, freversed.bond_ahead)
            if bond_cmp > 0:
                return straight, straight_features
            if bond_cmp < 0:
                return reverse, reverse_features

        return straight, straight_features

    def decompose_name(self, decomposition: MolDecomposition):
        features = list(self.features(decomposition))
