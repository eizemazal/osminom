import heapq
from collections import defaultdict, deque
from dataclasses import dataclass
from enum import Enum
from functools import cached_property, lru_cache
from typing import Generator, TypeAlias

from buftinom.funcgroups import GroupMatch, get_matchers
from buftinom.lookup import (
    BOND_PRIORITY,
)
from buftinom.smileg import Atom, Molecule
from buftinom.utils import filter_max

ChainKey: TypeAlias = tuple[Atom, Atom]
Chain: TypeAlias = list[Atom] | tuple[Atom, ...]


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
    functional_groups: dict[Atom, list[GroupMatch]]
    is_cycle: bool
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

    def with_chain(self, chain: Chain):
        return MolDecomposition(
            chain=chain,
            connections=self.connections,
            functional_groups=self.functional_groups,
            is_cycle=self.is_cycle,
            fpod=self.fpod,
        )


@dataclass
class ChainGroup:
    matcher: set[Atom]
    chains: list[Chain]


class Alogrythms:
    def __init__(self, mol: Molecule):
        self.mol = mol
        self.chains: list[Chain] = []

    @cached_property
    def leafs(self):
        if len(self.mol.atoms) == 1:
            return [self.mol.atoms[0]]

        cycles = self.cycles
        cycle_atoms = set()
        for cycle in cycles:
            cycle_atoms |= set(cycle)

        leafs = []
        for atom, connections in self.mol.adj.items():
            if atom in self.functional_group_atoms:
                continue
            if atom in cycle_atoms:
                continue

            available_conns = (
                set(connections)
                .difference(cycle_atoms)
                .difference(self.functional_group_atoms)
            )
            if len(available_conns) <= 1:
                leafs.append(atom)

        # Take one point from the cycle
        # so we can track chainded connections from the leafs and the cycles
        for cycle in cycles:
            leafs.append(cycle[0])

        return leafs

    @cached_property
    def functional_groups(self):
        matchers = get_matchers(self.mol)

        result: dict[Atom, GroupMatch] = {}

        for atom in self.mol.atoms:
            for matcher in matchers:
                mtch = matcher.matches(atom)
                if mtch and mtch.root not in result:
                    result[mtch.root] = mtch

        return result

    @cached_property
    def functional_group_atoms(self):
        atoms: set[Atom] = set()
        for func_group in self.functional_groups.values():
            group_atoms = func_group.atoms
            if func_group.root:
                group_atoms = group_atoms.difference({func_group.root})

            atoms |= group_atoms
        return atoms

    def cycle_exists(self, cycle_matchers: list[set[Atom]], cycle_atoms: set[Atom]):
        for matcher in cycle_matchers:
            if cycle_atoms == matcher:
                return True

        return False

    @cached_property
    def cycles(self) -> list[Chain]:
        """BFS/ modified Johnson's algorithm for undirected graph"""
        adj = self.mol.adj
        atoms = self.mol.atoms
        adjset = {atom: set(connections) for atom, connections in adj.items()}

        visited: set[Atom] = set()
        saturated: set[Atom] = set()

        cycles: list[Chain] = []
        cycle_matchers: list[set[Atom]] = []

        for atom in atoms:
            cycle = self.shortest_cycle(
                adj,
                atom,
                ignore_start=visited,
                always_ignore=saturated,
            )
            cycle_atoms = set(cycle)

            if not cycle:
                continue

            if self.cycle_exists(cycle_matchers, cycle_atoms):
                continue

            cycle_matchers.append(cycle_atoms)
            cycles.append(cycle)

            visited.add(atom)
            visited |= cycle_atoms

            for c in cycle + [atom]:
                if c in saturated:
                    continue

                if len(adjset[c].difference(visited)) == 0:
                    # No more visits for atoms which all connections
                    # Already been processed
                    saturated.add(c)

        return cycles

    def shortest_cycle(
        self,
        graph: dict[Atom, list[Atom]],
        start: Atom,
        ignore_start: set[Atom] = None,
        always_ignore: set[Atom] = None,
    ) -> Chain:
        """bfs"""
        if always_ignore is None:
            always_ignore = set()

        queue = deque([(start, None)])
        visited = {start}
        predecessor: dict[Atom, Atom | None] = {start: None}

        while queue:
            current_vertex, parent = queue.popleft()

            neighbors = graph[current_vertex]

            # do not start with paths that already been traversed
            if ignore_start:
                nbset = set(neighbors)
                allowed_starts = nbset.difference(ignore_start)
                neighbors = allowed_starts
                ignore_start = None

            for neighbor in neighbors:
                if neighbor in always_ignore:
                    continue

                if neighbor not in visited:
                    visited.add(neighbor)
                    #
                    predecessor[neighbor] = current_vertex
                    queue.append((neighbor, current_vertex))
                elif neighbor != parent and parent is not None:
                    # Found a cycle, now reconstruct it
                    return self.unfold_cycle(predecessor, neighbor, current_vertex)

        return []

    def unfold_cycle(self, predc: dict[Atom, Atom], start: Atom, end: Atom) -> Chain:
        """Since we use bfs cycle will look like a collar"""
        left = [end]
        right = [start]

        while predc[left[-1]] != None:
            left.append(predc[left[-1]])

        while predc[right[-1]] != None:
            right.append(predc[right[-1]])

        left = reversechain(left)
        right = reversechain(right)
        last_common_node = left[0]
        assert right[0] == last_common_node

        strip_tail_index = 0
        for idx, (l, r) in enumerate(zip(left, right)):
            if l == r:
                last_common_node = l
                continue

            strip_tail_index = idx
            break

        left = left[strip_tail_index:]
        right = right[strip_tail_index:]

        path = left + reversechain(right) + [last_common_node]
        return path

    def max_cycle(self, cycles: list[Chain]) -> Chain:
        maxlen = max(len(c) for c in cycles)
        by_length = [c for c in cycles if len(c) == maxlen]
        # if len(by_length) == 1:
        #     return by_length[0]
        # if aromatic - better
        # if func group - even better
        # now only len
        return by_length[0]

    @lru_cache(maxsize=256)
    def distances_from(self, start: Atom):
        """Dijkstra's algorithm"""
        adj = self.mol.adj
        queue = [(0, start)]

        distances = {a: float("inf") for a in self.mol.atoms}
        distances[start] = 0

        predecessors: dict[Atom, Atom] = {a: None for a in self.mol.atoms}

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

    def unfold_chain(self, predc: dict[Atom, Atom], start: Atom, end: Atom) -> Chain:
        path = []
        node = end
        while node is not None:
            path.append(node)
            node = predc[node]

        path.reverse()
        return path if path[0] == start else []

    @cached_property
    def leaf_distances(self):
        """Calculate distances between all leaf atoms"""
        leafs = self.leafs
        leaf_distances: dict[ChainKey, float] = {}

        if len(leafs) == 1:
            leaf = [leafs[0]]
            self.chains.append(leaf)
            return {chainkey(leaf): 0}

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

    def chain_count_functional_group(self, chain: Chain):
        i = 0
        for atom in chain:
            if atom in self.functional_groups:
                i += 1
        return i

    def chain_length(self, chain: Chain):
        return len(chain)

    def chain_bonds(self, chain: Chain):
        if len(chain) < 2:
            return 0

        prio = 0
        for a1, a2 in zip(chain, chain[1:]):
            bond = self.mol.bonds[(a1, a2)]
            prio += BOND_PRIORITY.get(bond.type, 0)
        return prio

    def chain_outside_connections(self, chain: Chain):
        chain_atoms = set(chain)
        nneighbors = 0
        for atom in chain:
            neighbors = set(self.mol.adj[atom]).difference(chain_atoms)
            nneighbors += len(neighbors)
        return nneighbors

    def max_chain(self, chains: list[Chain]) -> Chain:
        max_chains = filter_max(chains, self.chain_count_functional_group)
        max_chains = filter_max(max_chains, self.chain_length)
        max_chains = filter_max(max_chains, self.chain_bonds)
        max_chains = filter_max(max_chains, self.chain_outside_connections)
        #
        max_chain = max_chains[0]
        #
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
        """Split chains by in graph components they in"""
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

    def assert_not_implemented(self, cycles: list[Chain]):
        # level 1
        if len(cycles) > 1:
            raise AssertionError("Multiple cycles is not implemented")

        # level 2
        for i, c1 in enumerate(cycles):
            for j, c2 in enumerate(cycles):
                if i == j:
                    continue
                if set(c1) & set(c2):
                    raise AssertionError(
                        f"Intersecting cycles are not implemented \n{c1} \n{c2}"
                    )

    def chain_groups(self, chain: Chain):
        chain_atoms = set(chain)
        return {
            atom: group
            for atom, group in self.functional_groups.items()
            if atom in chain_atoms
        }

    def decompose(self):
        chains = self.all_chains
        cycles = self.cycles
        self.assert_not_implemented(cycles)

        def _decompose_cycles(chains: list[Chain], cycles: list[Chain]):
            if not cycles:
                return _decompose(chains)

            connections: dict[Atom, list[Chain]] = defaultdict(list)

            max_cycle = self.max_cycle(cycles)

            subchains = self.stripchains(chains, max_cycle)
            groupped = self.group_by_ends(subchains)

            for connection, groupped_chains in groupped.items():
                for chaingroup in self.interconneced(groupped_chains):
                    connections[connection].append(_decompose(chaingroup))

            groups = self.chain_groups(max_cycle)

            return MolDecomposition(
                chain=tuple(max_cycle),
                connections=connections,
                is_cycle=True,
                functional_groups=groups,
            )

        def _decompose(chains: list[Chain]):
            connections: dict[Atom, list[Chain]] = defaultdict(list)

            max_chain = self.max_chain(chains)
            subchains = self.stripchains(chains, max_chain)
            groupped = self.group_by_ends(subchains)

            for connection, groupped_chains in groupped.items():
                for chaingroup in self.interconneced(groupped_chains):
                    connections[connection].append(_decompose(chaingroup))

            groups = self.chain_groups(max_chain)

            return MolDecomposition(
                chain=tuple(max_chain),
                connections=connections,
                is_cycle=False,
                functional_groups=groups,
            )

        return _decompose_cycles(chains, cycles)
