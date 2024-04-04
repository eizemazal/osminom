import heapq
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property, lru_cache
from typing import TypeAlias

from buftinom.lookup import ROOT_BY_LENGTH, PrimarySuffix, bond_prio_cmp
from buftinom.smileg import Atom, BondType, Molecule, MoleculeConstructor

Chainmap: TypeAlias = dict[tuple[Atom, Atom], list[Atom]]


def chainkey(chain: list[Atom]) -> tuple[Atom, Atom]:
    return chain[0], chain[-1]


@dataclass
class MolDecomposition:
    chain: list[Atom]
    connections: dict[Atom, "MolDecomposition"]

    def print(self, level=0):
        print(self.chain)
        for connection, decomp in self.connections.items():
            print("    " * level, connection, ":", end=" ")
            decomp.print(level + 2)


class Alogrythms:
    def __init__(self, mol: Molecule):
        self.mol = mol
        self.chains: Chainmap = {}

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

        predecessors = {a: None for a in self.mol.atoms}

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

    def _unfold_chain(self, predecessors, start, end):
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
        leaf_distances: Chainmap = {}

        for a1 in leafs:
            a1_distances, a1_predecessors = self.distances_from(a1)

            for a2 in leafs:
                if a1 == a2:
                    continue
                self.chains[(a1, a2)] = self._unfold_chain(a1_predecessors, a1, a2)
                leaf_distances[(a1, a2)] = a1_distances[a2]

        return leaf_distances

    @cached_property
    def deduped_chains(self):
        self.leaf_distances

        deduped: Chainmap = {}

        for (a1, a2), chain in self.chains.items():
            if (a1, a2) in deduped or (a2, a1) in deduped:
                continue

            deduped[(a1, a2)] = chain

        return deduped

    def max_chains(self, chains: Chainmap) -> Chainmap:
        longest_chain_len = len(max(chains.values(), key=len))

        return {
            index: chain
            for index, chain in chains.items()
            if len(chain) == longest_chain_len
        }

    def max_chain(self, chains: Chainmap) -> list[Atom]:
        max_chains = self.max_chains(chains)
        assert len(max_chains) == 1, "Not implemented for multiple max chains"
        _, max_chain = max_chains.popitem()
        return max_chain

    def orient_by_leafs(self, chains: Chainmap) -> Chainmap:
        leafs = self.leafs
        oriented: Chainmap = {}

        for (a1, a2), chain in chains.items():
            if a2 in leafs:
                oriented[(a2, a1)] = list(reversed(chain))
            else:
                oriented[(a1, a2)] = chain

        return oriented

    def stripchains(self, chains: Chainmap, chain: list[Atom]):
        """Generate all subchains
        All subchains will start with leaf, and end on the input chain
        """
        chainset = set(chain)
        chains = self.orient_by_leafs(chains)
        splits: Chainmap = {}

        for splitting in chains.values():
            buffer: list[Atom] = []
            for a in splitting:
                if a in chainset:
                    if buffer:
                        buffer.append(a)
                        splits[chainkey(buffer)] = list(buffer)
                    buffer = []
                    continue

                buffer.append(a)

        return splits

    def group_by_ends(self, chains: Chainmap):
        grouped: dict[Atom, Chainmap] = defaultdict(dict)

        for (start, end), chain_end in chains.items():
            assert len(chain_end) > 1, chain_end
            *chain, end = chain_end

            # reverse to orient the chain from the groupping atom forward
            grouped[end][(chain[-1], start)] = list(reversed(chain))

        return grouped

    def decompose(self):
        chains = self.deduped_chains

        def _decompose(chains):
            max_chain = self.max_chain(chains)
            subchains = self.stripchains(chains, max_chain)
            groupped = self.group_by_ends(subchains)

            connections = {}
            for connection, chains in groupped.items():
                connections[connection] = _decompose(chains)

            return MolDecomposition(max_chain, connections)

        return _decompose(chains)


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
            return (0, (self.chain[0], BondType.SINGLE, self.chain[0]))

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


class Iupac:
    def __init__(self, mol: MoleculeConstructor):
        self.mol = mol

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
                form = suffix.value.norm

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
