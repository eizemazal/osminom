from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property, partial
from typing import Callable, NamedTuple

from buftinom.algorythms import (
    Alogrythms,
    Chain,
    MolDecomposition,
    reversechain,
)
from buftinom.funcgroups import GroupMatch
from buftinom.lookup import (
    MULTI_BY_PREFIX,
    MULTI_MULTI_BY_PREFIX,
    ROOT_BY_LENGTH,
    Alphabet,
    Infix,
    Prefix,
    PrimarySuffix,
)
from buftinom.smileg import Atom, Bond, BondType, Molecule
from buftinom.translate import WordForm
from buftinom.utils import choose


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


def singleidx2str(name: WordForm, idx: int):
    if idx > 1:
        return ["-", str(idx), "-", name.short]
    else:
        return [name.short]


def manyidx2str(name: WordForm, ids: list[int]):
    index = ",".join(map(str, ids))
    multi = MULTI_BY_PREFIX.get(len(ids))
    return ["-", index, "-", multi.value.norm, name.short]


def suffixes2str(suffixes: list[Synt], sub_suffix: WordForm):
    res: list[str] = []

    ssuff: list[Synt] = sorted(suffixes, key=lambda s: s.value)

    groupped: dict[str, list[WordForm]] = defaultdict(list)
    for idx, name in suffixes:
        groupped[name].append(idx)

    groups = sorted(groupped.items(), key=lambda g: g[0])
    for name, ids in groups:
        if len(ids) > 1:
            res.extend(manyidx2str(name, ids))
        else:
            res.extend(singleidx2str(name, ids[0]))

    if sub_suffix:
        res.extend(sub_suffix.norm)
    else:
        res[-1] = ssuff[-1].value.norm

    return "".join(res)


def preffixes2str(iupac: IupacName):
    preffixes = iupac.prefixes
    if not preffixes:
        return None

    #
    subnames = [(idx, iupac2str(subn)) for idx, subn in preffixes]
    #

    groups: dict[str, list[int]] = defaultdict(list)

    for idx, name in subnames:
        groups[name].append(idx)

    ordered = sorted(list(groups.items()), key=lambda g: g[0])

    res: list[str] = []

    for name, indexes in ordered:
        subname: list[str] = []

        name_complex = name[0].isdigit()
        multiselector = MULTI_MULTI_BY_PREFIX if name_complex else MULTI_BY_PREFIX

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


@dataclass
class AtomFeatures:
    chain_index: int
    atom: Atom
    bond_ahead: Bond
    subciupacs: list[IupacName]
    functional_group: GroupMatch | None


class Iupac:
    def __init__(self, mol: Molecule):
        self.mol = mol

    @cached_property
    def decomposition(self):
        return Alogrythms(self.mol).decompose()

    def features(
        self, decomp: MolDecomposition, subiupacs: dict[Atom, list[IupacName]]
    ):
        chain = decomp.chain

        if decomp.is_cycle:
            chain: Chain = decomp.chain + (decomp.chain[0],)

        for i, (a1, a2) in enumerate(zip(chain, chain[1:]), start=1):
            bond = self.mol.bonds[(a1, a2)]

            yield AtomFeatures(
                chain_index=i,
                atom=a1,
                bond_ahead=bond,
                subciupacs=subiupacs.get(a1),
                functional_group=decomp.functional_groups.get(a1),
            )

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

    def preffixes(self, dec: MolDecomposition):
        result: dict[Atom, list[IupacName]] = defaultdict(list)
        for atom in dec.chain:
            subchains = dec.connections.get(atom, [])

            for subchain in subchains:
                #
                subiupac = self.decompose_name(subchain, primary=False)
                #
                result[atom].append(subiupac)

        return result

    def index_preffixes(
        self, features: list[AtomFeatures], all_preffixes: dict[Atom, list[IupacName]]
    ):
        result = []

        for feature in features:
            if preffixes := all_preffixes.get(feature.atom):
                for preffix in preffixes:
                    result.append(Subn(feature.chain_index, preffix))
        return result

    PRIORITIES = {
        BondType.SINGLE: 100,
        "weak-functional-group": 200,
        BondType.TRIPLE: 300,
        BondType.DOUBLE: 400,
        "functional-group": 1000,
    }

    def feature_bonds(self, feature: AtomFeatures):
        return Iupac.PRIORITIES.get(feature.bond_ahead.type)

    def feature_group(self, feature: AtomFeatures):
        prio = Iupac.PRIORITIES
        if feature.functional_group is not None:
            group_tag = feature.functional_group
            if group_tag in {}:
                return 0
            return prio.get("functional-group")
        return 0

    def feature_weak_group(self, feature: AtomFeatures):
        prio = Iupac.PRIORITIES
        if feature.functional_group is not None:
            group_tag = feature.functional_group
            if group_tag in {}:
                return prio.get("weak-functional-group")

        return 0

    def feature_subchain(self, feature: AtomFeatures):
        alphabetical = 0
        if feature.subciupacs:
            for sub in feature.subciupacs:
                # handle for complex radicals
                alphabetical = max(alphabetical, ord(sub.root_word.norm[0]))

            alphabetical = -alphabetical + ord(Alphabet.MAX.value.norm)
        return alphabetical

    def compare(
        self,
        *,
        fs1: list[AtomFeatures],
        fs2: list[AtomFeatures],
        key: Callable[[AtomFeatures], float],
    ):
        mfs1 = map(key, fs1)
        mfs2 = map(key, fs2)
        for score1, score2 in zip(mfs1, mfs2):
            if score1 != score2:
                return score1 - score2

        return 0

    def fpod(self, decomp: MolDecomposition, preffixes: dict[Atom, list[IupacName]]):
        straight = decomp
        reversed = decomp.with_chain(tuple(reversechain(straight.chain)))

        straightf = list(self.features(straight, preffixes))
        reversef = list(self.features(reversed, preffixes))

        selector = partial(choose, a=(straight, straightf), b=(reversed, reversef))
        comparator = partial(self.compare, fs1=straightf, fs2=reversef)

        if result := selector(cmp=comparator(key=self.feature_group)):
            return result

        if result := selector(cmp=comparator(key=self.feature_bonds)):
            return result

        if result := selector(cmp=comparator(key=self.feature_weak_group)):
            return result

        if result := selector(cmp=comparator(key=self.feature_subchain)):
            return result

        return straight, straightf

    def decompose_name(self, decomp: MolDecomposition, *, primary: bool = True):
        #
        unordered_preffixes = self.preffixes(decomp)
        #
        decomp, features = self.fpod(decomp, unordered_preffixes)

        preffixes = self.index_preffixes(features, unordered_preffixes)
        root = ROOT_BY_LENGTH[len(decomp.chain)].value
        suffix = self.suffixes_by_features(features, primary=primary)
        infix = Infix.CYCLO.value if decomp.is_cycle else None
        subsuff = Prefix.YL.value if not primary else None

        return IupacName(
            prefixes=preffixes,
            infix=infix,
            root_word=root,
            prime_suffixes=suffix,
            sub_suffix=subsuff,
            func_suffix=None,
            ref=decomp,
        )
