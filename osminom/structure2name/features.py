from collections import defaultdict
from dataclasses import dataclass, field

from osminom.structure2name.algorithms import (
    Chain,
    MolDecomposition,
    SubchainConnection,
    reversechain,
)
from osminom.structure2name.funcgroups import GroupMatch
from osminom.structure2name.lookup import Alphabet, FunctionalGroup
from osminom.structure2name.models import IupacName
from osminom.structure2name.molecule import Atom, Bond, BondType
from osminom.structure2name.utils import first_max, nonzero_indexes


@dataclass
class AtomFeatures:
    chain_index: int
    atom: Atom
    bond_ahead: Bond
    # self atom in this connection is peer
    connection: SubchainConnection
    subiupacs: list[IupacName]
    functional_groups: list[GroupMatch]


def f2c(features: list[AtomFeatures]) -> Chain:
    return tuple(f.atom for f in features)


@dataclass
class DecFeatures:
    dec: MolDecomposition
    fs: list[AtomFeatures]

    scores: dict[str, list[int]] = field(default_factory=lambda: defaultdict(list))

    def as_tuple(self):
        return self.dec, self.fs


def feature_connected(feature: AtomFeatures):
    return int(feature.connection is not None)


def feature_noncarbon(feature: AtomFeatures):
    return int(feature.atom.symbol.lower() != "c")


def feature_bonds(feature: AtomFeatures):
    if feature.bond_ahead is None:
        return 0
    return Features.PRIORITIES.get(feature.bond_ahead.type, 0)


def feature_group(feature: AtomFeatures):
    prio = Features.PRIORITIES
    if not feature.functional_groups:
        return 0

    score = 0

    for group in feature.functional_groups:
        group_tag = group.tag
        if group_tag not in Features.WEAK_FUNCTIONAL_GROUPS:
            score = prio.get("functional-group")

    return score


def feature_weak_group(feature: AtomFeatures):
    prio = Features.PRIORITIES
    if not feature.functional_groups:
        return 0

    score = 0
    for group in feature.functional_groups:
        group_tag = group.tag
        if group_tag in Features.WEAK_FUNCTIONAL_GROUPS:
            score = prio.get("functional-group")

    return score


def feature_subchain(feature: AtomFeatures):
    alphabetical = 0
    if feature.subiupacs:
        for sub in feature.subiupacs:
            # handle for complex radicals
            alphabetical = max(alphabetical, ord(sub.root_word.norm[0]))

        alphabetical = -alphabetical + ord(Alphabet.MAX.value.norm)
    return alphabetical


def calc_feature(feature_func, decs_features: list[DecFeatures]):
    for decfs in decs_features:
        decfs.scores[feature_func.__name__] = list(map(feature_func, decfs.fs))

    i = first_max([d.scores[feature_func.__name__] for d in decs_features])
    if i is not None:
        return decs_features[i].as_tuple()


class Features:
    PRIORITIES = {
        BondType.SINGLE: 0,
        "weak-functional-group": 1,
        BondType.TRIPLE: 2,
        BondType.DOUBLE: 3,
        "functional-group": 4,
    }

    WEAK_FUNCTIONAL_GROUPS = {
        FunctionalGroup.OXY,
        FunctionalGroup.AMINO,
        FunctionalGroup.ESTER,
        FunctionalGroup.BROM,
        FunctionalGroup.CHLOR,
        FunctionalGroup.NITRO,
    }

    def __init__(self, mol: MolDecomposition):
        self.mol = mol

    def features(
        self, decomp: MolDecomposition, subiupacs: dict[Atom, list[IupacName]]
    ):
        """
        Collect features of the decomposed chain.

        Spectial attention to cycles, their representation contains bonds between the edges of the chain
        """
        chain = decomp.chain
        parent_connection = decomp.connected_by

        if decomp.is_cycle:
            chain: Chain = decomp.chain + (decomp.chain[0],)

        non_carbons = 0
        for i, (a1, a2) in enumerate(zip(chain, chain[1:]), start=1):
            if a1.symbol.lower() != "c" and not decomp.is_aromatic:
                non_carbons += 1
                continue

            bond = self.mol.bonds[(a1, a2)]

            parent = None
            if parent_connection and a1 == parent_connection.peer:
                parent = decomp.connected_by

            func_groups = decomp.functional_groups.get(a1, [])
            # if not is_carbon(a2):
            #     nc_group = decomp.functional_groups.get(a2)
            #     if func_group is not None and nc_group is not None:
            #         raise NotImplementedError(f"Multiple func groups on one atom {a1}")

            #     func_group = nc_group

            index = i - non_carbons
            yield AtomFeatures(
                chain_index=index,
                atom=a1,
                bond_ahead=bond,
                subiupacs=subiupacs.get(a1),
                functional_groups=func_groups,
                connection=parent,
            )

        # Do not forget the last atom of the chain
        if decomp.is_cycle:
            return

        last = chain[-1]

        parent = None
        if parent_connection and last == parent_connection.peer:
            parent = decomp.connected_by

        yield AtomFeatures(
            chain_index=len(chain),
            atom=last,
            bond_ahead=None,
            subiupacs=subiupacs.get(last),
            functional_groups=decomp.functional_groups.get(last, []),
            connection=parent,
        )

    def fpod(self, decomp: MolDecomposition, preffixes: dict[Atom, list[IupacName]]):
        """
        First Point of Difference

        The algorithm to determine the Numbering of the chain.

        Short rule is:
        1. Primary functional group is the king, it have the least number
        2. Double and Triple bonds have second prio
        3. Weak functional groups next
        4. Afther that - Subchains.
        5. Other collisions resolved alphabetically

        We test every possible chain numberings to see which higher priority is matches
        And select numbering appropriately.

        For the acyclic chain all possibilities are - straight or reversed numbering
        """
        if decomp.is_cycle:
            return self.fpod_cycle(decomp, preffixes)

        straight = decomp
        reversed = decomp.with_chain(tuple(reversechain(straight.chain)))

        return self._fpod_many(preffixes, [straight, reversed])

    def fpod_cycle(
        self, decomp: MolDecomposition, preffixes: dict[Atom, list[IupacName]]
    ):
        """
        Cycles is a special case.
        We first look for any priority features,

        The fpod then is preformed with this feature as the first atom of the chain.
        And the position wins if in the end it'll have the best priority fit.
        """
        features = list(self.features(decomp, preffixes))

        priorities = [
            feature_connected,
            feature_noncarbon,
            feature_group,
            feature_bonds,
            feature_weak_group,
            feature_subchain,
            # any
            lambda _: 1,
        ]

        starts = []
        for key_func in priorities:
            starts = list(nonzero_indexes(features, key_func))
            if starts:
                break

        comparable = []

        #                   we are here
        # < here for backward  |  here for forward >
        #
        # test all possible starts in circle
        chain = decomp.chain + decomp.chain + decomp.chain

        for first_id in starts:
            curr_idx = len(decomp.chain)
            start = curr_idx + first_id
            size = len(decomp.chain)

            fchain = chain[start : start + size]
            bchain = reversechain(chain[start - size + 1 : start + 1])

            assert len(fchain) == len(bchain)

            comparable.append(decomp.with_chain(fchain))
            comparable.append(decomp.with_chain(bchain))

        return self._fpod_many(preffixes, comparable)

    def _fpod_many(
        self,
        preffixes: dict[Atom, list[IupacName]],
        decs: list[MolDecomposition],
    ):
        """
        For each given decomposition, compare features that it provides
        to select decomposition with the **first point of dfference**

        Example: assume we have CCC=C=C
        Considering two names pent-1,2-diene / pent-3,4-diene
        Bond scores will look like:
        [0, 0, 2, 2] - for CC=C=C
        [0, 2, 2, 0] - for reversed C=C=CC

        algo will go through the columns of this table and
        find the first column that will signal that it have score with highest prio

        in the examle it'll find that the best fit is second column second row,
        and we'll select reversed decomposition.
        """
        features = [
            DecFeatures(dec, list(self.features(dec, preffixes))) for dec in decs
        ]

        funcs = [
            feature_connected,
            feature_group,
            feature_bonds,
            feature_weak_group,
            feature_subchain,
        ]

        for feature_func in funcs:
            if selected_prio := calc_feature(feature_func, features):
                return selected_prio

        return features[0].as_tuple()
