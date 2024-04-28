import heapq
from collections import defaultdict, deque
from dataclasses import dataclass
from functools import cached_property, lru_cache, partial
from typing import Generator, TypeAlias

from buftinom.funcgroups import GroupMatch, get_matchers
from buftinom.lookup import (
    BOND_PRIORITY,
)
from buftinom.smileg import Atom, BondType, Molecule, is_carbon
from buftinom.utils import deepmerge, filter_max

ChainKey: TypeAlias = tuple[Atom, Atom]
Chain: TypeAlias = tuple[Atom, ...]


def chainkey(chain: Chain) -> ChainKey:
    assert len(chain) > 0, "Empty chain has no key"
    return chain[0], chain[-1]


def reversekey(key: ChainKey) -> ChainKey:
    return key[1], key[0]


def reversechain(chain: Chain) -> Chain:
    return tuple(reversed(chain))


@dataclass(eq=True, frozen=True)
class SubchainConnection:
    parent: Atom
    peer: Atom
    bond: BondType
    via: GroupMatch | None

    def __hash__(self):
        return hash((self.parent, self.peer, self.bond))

    def __str__(self):
        via = self.via
        if not via:
            via = self.bond
        return f"<{self.parent}{via}{self.peer}>"

    __repr__ = __str__


@dataclass
class MolDecomposition:
    chain: Chain
    connections: dict[SubchainConnection, list["MolDecomposition"]]
    functional_groups: dict[Atom, list[GroupMatch]]
    # if it is subchain we indicate to which atom it is connected
    # Here, the peer will be the atom from this chain
    connected_by: SubchainConnection | None
    is_cycle: bool
    is_aromatic: bool

    def extra_tags(self):
        tags = []
        if self.is_cycle:
            tags.append("cycle")
        if self.is_aromatic:
            tags.append("aromatic")
        return ",".join(tags)

    def print(self, level=0):
        print(self.chain, self.extra_tags())
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
            chain=tuple(chain),
            connections=self.connections,
            functional_groups=self.functional_groups,
            is_cycle=self.is_cycle,
            is_aromatic=self.is_aromatic,
            connected_by=self.connected_by,
        )


class Alogrythms:
    """
    A class that represents algorithms for working with molecules.

    This class provides various algorithms for analyzing and manipulating molecules.
    It includes methods for finding cycles, calculating distances, identifying functional groups,
    and generating subchains, among others.
    """

    def __init__(self, mol: Molecule):
        self.mol = mol
        self.chains: list[Chain] = []

    @cached_property
    def leafs(self):
        """Find 'leafs' of the molecule.

        An atom considered a leaf if it is
        - end of the chain,
        - not the atom of the functional group
        - atom of the cycle (only one needed)

        This leafs will be the bounds of the primary and side chains
        """
        if len(self.mol.atoms) == 1:
            return [self.mol.atoms[0]]

        cycle_atoms = self.cycle_atoms

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

        return leafs

    def on_cycle_border(self, atom: Atom):
        if atom in self.cycle_atoms:
            return {atom}

        if not atom:
            return set()
        return self.mol.adj_set[atom] & self.cycle_atoms

    @cached_property
    def functional_groups(self):
        """Matching and finding functional groups of the molecule

        Matching is ordered, so the first found group of the atom
        will take priority

        Limitations: one functional group to one atom
        """
        matchers = get_matchers(self.mol)

        groups: dict[Atom, list[GroupMatch]] = defaultdict(list)
        group_atoms: set[Atom] = set()

        for matcher in matchers:
            for atom in self.mol.atoms:
                mtch = matcher.matches(atom)
                if not mtch:
                    continue

                match_atom_used = mtch.func_atoms & group_atoms

                if match_atom_used:
                    continue

                if is_carbon(mtch.root) and mtch.is_flex_root:
                    if croot := self.on_cycle_border(mtch.root):
                        if len(croot) != 1:
                            raise NotImplementedError(
                                f"Atom {mtch.root} have too much connections with the cycle {croot}"
                            )
                        group_atoms |= {mtch.root}
                        mtch.change_root(croot.pop())

                groups[mtch.root].append(mtch)
                if not mtch.is_symmetric:
                    group_atoms |= mtch.func_atoms

        return groups

    @cached_property
    def functional_group_atoms(self):
        """
        Lil helper, returns atoms of the functional groups, excluding carbon atoms

        This atoms used to ignore them during pathfinding
        """
        atoms: set[Atom] = set()
        for func_groups in self.functional_groups.values():
            for group in func_groups:
                atoms |= group.func_atoms

                if group.root_flexed:
                    atoms |= {group.root_flexed}

        nonfunc_atoms = set(self.mol.atoms) - atoms
        noncarbon_func = list(filter(lambda a: not is_carbon(a), nonfunc_atoms))
        if noncarbon_func:
            raise NotImplementedError(
                f"Atoms {noncarbon_func} is possibly a part of the"
                + " functional group that is not yet supported"
            )

        return atoms

    def cycle_exists(self, cycle_matchers: list[set[Atom]], cycle_atoms: set[Atom]):
        for matcher in cycle_matchers:
            if cycle_atoms == matcher:
                return True

        return False

    @cached_property
    def cycles(self) -> list[Chain]:
        """BFS/ modified Johnson's algorithm for undirected graph

        Find **all** cycles of the atom. It will be minimal-sized cycles
        without nested cycles (cycles that contain other cycles inside)
        """
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

            for c in cycle + (atom,):
                if c in saturated:
                    continue

                if len(adjset[c].difference(visited)) == 0:
                    # No more visits for atoms which all connections
                    # Already been processed
                    saturated.add(c)

        return cycles

    @cached_property
    def cycle_atoms(self):
        cycle_atoms: set[Atom] = set()
        for cycle in self.cycles:
            cycle_atoms |= set(cycle)

        return cycle_atoms

    def shortest_cycle(
        self,
        graph: dict[Atom, list[Atom]],
        start: Atom,
        ignore_start: set[Atom] = None,
        always_ignore: set[Atom] = None,
    ) -> Chain:
        """bfs, shortest cycle

        Takes additional params, used by self.cycles to modify search to ignore already found cycles
        """
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
        """Since we use bfs cycle will look like a collar, make it a ring"""
        left = [end]
        right = [start]

        while predc[left[-1]] is not None:
            left.append(predc[left[-1]])

        while predc[right[-1]] is not None:
            right.append(predc[right[-1]])

        left = reversechain(left)
        right = reversechain(right)
        last_common_node = left[0]
        assert right[0] == last_common_node

        strip_tail_index = 0
        for idx, (latom, ratom) in enumerate(zip(left, right)):
            if latom == ratom:
                last_common_node = latom
                continue

            strip_tail_index = idx
            break

        left = left[strip_tail_index:]
        right = right[strip_tail_index:]

        path = left + reversechain(right) + (last_common_node,)
        return path

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

    def unfold_chain(
        self, predc: dict[Atom, Atom], start: Atom, end: Atom, ignore: set[Atom]
    ) -> Chain:
        path = []
        node = end
        while node is not None:
            if node in ignore:
                break
            path.append(node)
            node = predc[node]

        path.reverse()
        return tuple(path)

    @cached_property
    def leaf_distances(self):
        """Calculate distances between all leaf atoms, using self.distances_from"""
        leafs = self.leafs
        leaf_distances: dict[ChainKey, float] = {}
        path_stops = self.cycle_atoms | self.functional_group_atoms

        if len(leafs) == 1:
            leaf = [leafs[0]]
            self.chains.append(leaf)
            return {chainkey(leaf): 0}

        for a1 in leafs:
            a1_distances, a1_predecessors = self.distances_from(a1)

            for a2 in leafs:
                if a1 == a2:
                    continue
                chain = self.unfold_chain(a1_predecessors, a1, a2, ignore=path_stops)
                if chain:
                    self.chains.append(chain)
                    leaf_distances[(a1, a2)] = a1_distances[a2]

        return leaf_distances

    @cached_property
    def all_chains(self) -> list[Chain]:
        self.leaf_distances
        if len(self.mol.atoms) == 1:
            a = self.mol.atoms[0]
            return [(a,)]

        return self.chains

    ## Begin score functions to calculate the priorities for selecting primary chain

    def chain_have_connection(
        self, connection: SubchainConnection | None, chain: Chain
    ):
        if connection is None:
            return 0

        atom_conn = connection.peer
        for atom in chain:
            if atom == atom_conn:
                return 1
        return 0

    def chain_count_functional_group(self, chain: Chain):
        i = 0
        for atom in chain:
            if atom in self.functional_groups:
                # one is enough :)
                i = 1
        return i

    def is_aromatic(self, cycle: Chain):
        """
        given the cycle chain (n > 2, c[0] == c[-1]) returns wether this cycle is aromatic
        """

        # Case 1: Hueckel's 4N+2 criterion by lowercase notation.
        # Raises error due if it is not true - it is a syntax error
        num_aroma = sum([a.symbol.islower() for a in cycle])
        if num_aroma > 0:
            if num_aroma > 2 and (num_aroma - 2) % 4 == 0:
                return True

            raise ValueError(
                "Invalid aromatic bonds specification "
                + "(4N+2 atoms should be in aromatic cycle)"
                + f"\nGot: {cycle}"
            )

        # Filter Hueckel's 4N+2 criteria. Test if cycle could possibly be aromatic
        n = len(cycle)
        if (n - 2) % 4 != 0:
            return False

        closure = cycle + (cycle[0],)

        # Case 2: all the atoms connected explicitly via `:`
        all_aromatic = True
        for a1, a2 in zip(closure, closure[1:]):
            if self.mol.bonds[(a1, a2)].type != BondType.AROMATIC:
                all_aromatic = False
                break

        if all_aromatic:
            return True

        # Case 3: Traditional C-C=C-C=C-... aromatic cycle definition
        expected = {BondType.SINGLE, BondType.DOUBLE}
        next_bond = {
            BondType.SINGLE: {BondType.DOUBLE},
            BondType.DOUBLE: {BondType.SINGLE},
        }
        mix_correct = True
        for a1, a2 in zip(closure, closure[1:]):
            bond = self.mol.bonds[a1, a2]
            if bond.type not in expected:
                mix_correct = False
                break

            expected = next_bond.get(bond.type)

        return mix_correct

    def chain_cycle(self, chain: Chain):
        if len(chain) < 3:
            return False
        first, last = chainkey(chain)
        return (last, first) in self.mol.bonds

    def chain_cycle_aromatic(self, chain: Chain):
        return self.chain_cycle(chain) and self.is_aromatic(chain)

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

    ## End score functions

    def max_chain(self, chains: list[Chain], connection: Atom = None) -> Chain:
        """
        Scores the chains by the priorities, and filter untill the main chain is selected
        """
        max_chains = chains
        filters = [
            partial(self.chain_have_connection, connection),
            self.chain_count_functional_group,
            self.chain_cycle,
            self.chain_length,
            self.chain_cycle_aromatic,
            self.chain_bonds,
            self.chain_outside_connections,
        ]

        for filter in filters:
            max_chains = filter_max(max_chains, filter)

        #
        max_chain = max_chains[0]
        #
        return max_chain

    def orient_by_leafs(self, chains: list[Chain]) -> list[Chain]:
        """
        Rotate chains, so the leaf atom will be the first atom

        Usefull to find connections with other chains
        """
        leafs = self.leafs
        oriented: list[Chain] = []

        for chain in chains:
            end = chain[-1]
            if end in leafs:
                chain = reversechain(chain)

            oriented.append(chain)

        return oriented

    def stripchains(self, chains: list[Chain], chain: Chain):
        """Generate all subchains.
        All subchains will start with leaf, and end on the input chain

        Removes atoms of param:chain from the param:chains,
        but keeps the atom that connects chains
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
                break

            if not buffer:
                continue

            if self.chain_cycle(buffer):
                splits.append(tuple(buffer))
                continue

            key = chainkey(buffer)
            if key not in existing:
                splits.append(tuple(buffer))
                existing.add(key)
                existing.add(reversekey(key))

        return splits

    def connection_point(self, group: list[Chain], base_chain: Chain):
        """
        Reutrn Atom by which this group connected to the chain
        """
        adj = self.mol.adj_set
        chain_atoms = set(base_chain)

        # base_chain to group atom
        connections: set[(Atom, Atom)] = set()

        for chain in group:
            for atom in chain:
                for e in adj[atom] & chain_atoms:
                    connections.add((e, atom))
                    break

        if len(connections) > 1:
            raise NotImplementedError(f"Bicyclo case, {connections}")

        if len(connections) == 1:
            base_atom, group_atom = connections.pop()

            return SubchainConnection(
                parent=base_atom,
                peer=group_atom,
                bond=self.mol.bonds[(base_atom, group_atom)].type,
                via=None,
            )

        # Check connection through the functional group (i.e. ester R-COO-R)

        # group match connections
        connections: set[GroupMatch] = set()

        for atom, func_groups in self.chain_functional_groups(base_chain).items():
            for func_group in func_groups:
                if func_group.side_root is None:
                    continue

                for chain in group:
                    for a in chain:
                        if func_group.side_root == a:
                            connections.add(func_group)

        if len(connections) > 1:
            raise NotImplementedError(f"Bicyclo case, {connections}")

        if len(connections) == 1:
            match = connections.pop()

            return SubchainConnection(
                parent=match.root,
                peer=match.side_root,
                bond=BondType.SINGLE,
                via=match,
            )

        # No direct connections
        return None

    def group_connections(self, groups: list[list[Chain]], chain: Chain):
        """
        for each group return atom by which this group connected to the chain
        """
        res: list[tuple[SubchainConnection, list[Chain]]] = []
        connectles_group: list[Chain] = []
        for group in groups:
            connector = self.connection_point(group, chain)
            if connector:
                res.append((connector, group))
                continue

            connectles_group.extend(group)

        return res, connectles_group

    def chain_in_component(self, chain: Chain, component: set[Atom]):
        for atom in chain:
            if atom in component:
                return True

    def group_by_components(self, chains: list[Chain], components: list[set[Atom]]):
        """
        Return chains components by graph components formed by stripping chain from the graph
        """
        groups: list[list[Chain]] = []
        for _ in range(len(components)):
            groups.append([])

        for chain in chains:
            for i, component in enumerate(components):
                if self.chain_in_component(chain, component):
                    groups[i].append(chain)
                    break

        return [g for g in groups if g]

    def components(self, chains: list[Chain]):
        """
        Return components of the graph formed by stripping chain
        """
        atoms: set[Atom] = set()
        for chain in chains:
            atoms |= set(chain)

        visited: set[Atom] = set()
        adj = self.mol.adj

        def dfs(atom: Atom):
            visited.add(atom)
            component = {atom}

            for neighbor in adj[atom]:
                if neighbor in atoms and neighbor not in visited:
                    component |= dfs(neighbor)

            return component

        components = []
        for atom in atoms:
            if atom not in visited:
                components.append(dfs(atom))

        return components

    def assert_not_implemented(self, cycles: list[Chain]):
        # level 1
        # if len(cycles) > 1:
        #     raise AssertionError("Multiple cycles are not implemented")

        # level 2
        for i, c1 in enumerate(cycles):
            for j, c2 in enumerate(cycles):
                if i == j:
                    continue
                if set(c1) & set(c2):
                    raise NotImplementedError(
                        f"Intersecting cycles are not implemented \n{c1} \n{c2}"
                    )

    def chain_functional_groups(self, chain: Chain):
        """
        Map the found functional groups to the root atom they are connected to
        """
        chain_atoms = set(chain)
        return {
            atom: group
            for atom, group in self.functional_groups.items()
            if atom in chain_atoms
        }

    def flex_groups(self, decomp: MolDecomposition, connection: SubchainConnection):
        ffg: dict[Atom, list[GroupMatch]] = defaultdict(list)

        if len(decomp.chain) != 1:
            return ffg

        for fg in decomp.functional_groups.get(decomp.chain[0], []):
            if fg.is_flex_root:
                fg.change_root(connection.parent)
                ffg[connection.parent].append(fg)
                continue

            ffg = {}
            break

        return ffg

    def all_atoms(self, decomp: MolDecomposition):
        atoms = set(decomp.chain)
        for groups in decomp.functional_groups.values():
            for group in groups:
                atoms |= group.atoms

        for conns in decomp.connections.values():
            for conn in conns:
                atoms |= self.all_atoms(conn)

        return atoms

    def validate(self, decomp: MolDecomposition):
        atoms = self.all_atoms(decomp)
        mol_atoms = set(self.mol.atoms)

        if atoms != mol_atoms:
            raise ValueError(
                "Cannot figure molecule structure. "
                + f"Unable to connect {mol_atoms - atoms} to the chain."
                + "\nMolecule, is invalid or, more likely, not supported 🫡"
            )

    def decompose(self):
        """
        Crate decomposition of the Molecule

        This is the main payload of the class.
        Creates MoleculeDecomposition object that contains recursive
        information about Molecule connections, functional groups and subchains.

        Decomposition selects IUPAC-correct primary chain and subchains.
        """
        chains = self.all_chains
        cycles = self.cycles
        self.assert_not_implemented(cycles)

        def _decompose(chains: list[Chain], parent_connection: SubchainConnection):
            connections: dict[Atom, list[Chain]] = defaultdict(list)

            max_chain = self.max_chain(chains, parent_connection)
            subchains = self.stripchains(chains, max_chain)
            components = self.components(subchains)
            chain_groups = self.group_by_components(subchains, components)
            groups, connectles = self.group_connections(chain_groups, max_chain)

            flex_func_groups: dict[Atom, list[GroupMatch]] = defaultdict(list)

            for connection, groupped_chains in groups:
                # Add all the chains without connections, to test
                # if they are connected to current subgroup
                summol_chains = groupped_chains + connectles
                subdecomp = _decompose(summol_chains, connection)

                if ffg := self.flex_groups(subdecomp, connection):
                    flex_func_groups = deepmerge(flex_func_groups, ffg)
                    continue

                connections[connection].append(subdecomp)

            func_groups = self.chain_functional_groups(max_chain)
            func_groups = deepmerge(func_groups, flex_func_groups)

            is_cycle = self.chain_cycle(max_chain)
            is_aromatic = is_cycle and self.is_aromatic(max_chain)

            return MolDecomposition(
                chain=tuple(max_chain),
                connections=connections,
                is_cycle=is_cycle,
                is_aromatic=is_aromatic,
                connected_by=parent_connection,
                functional_groups=func_groups,
            )

        final_decomposition = _decompose(cycles + chains, None)
        self.validate(final_decomposition)

        return final_decomposition
