from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from typing import Literal, NamedTuple

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
from buftinom.utils import first_max, nonzero_indexes


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
    func_suffixes: list[Synt]

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


def singleidx2str(name: WordForm, idx: int, form: Literal["short", "norm"]):
    if idx > 1:
        return ["-", str(idx), "-", name.get(form)]
    else:
        return [name.get(form)]


def manyidx2str(name: WordForm, ids: list[int], form: Literal["short", "norm"]):
    index = ",".join(map(str, ids))
    multi = MULTI_BY_PREFIX.get(len(ids))
    return ["-", index, "-", multi.value.norm, name.get(form)]


def suffixes2str(suffixes: list[Synt], form: Literal["short", "norm"]):
    if not suffixes:
        return [], None

    res: list[str] = []

    groupped: dict[WordForm, list[int]] = defaultdict(list)
    for idx, name in suffixes:
        groupped[name].append(idx)

    groups = sorted(groupped.items(), key=lambda g: g[0])
    lastname = None
    for name, ids in groups:
        if len(ids) > 1:
            res.extend(manyidx2str(name, ids, form))
        else:
            res.extend(singleidx2str(name, ids[0], form))

        lastname = name

    return res, lastname


def all_suffixes2str(iupac: IupacName):
    primes, lastprime = suffixes2str(iupac.prime_suffixes, form="short")
    functionals, _ = suffixes2str(iupac.func_suffixes, form="norm")

    res = primes

    if iupac.sub_suffix:
        res.append(iupac.sub_suffix.norm)
    if functionals:
        res.extend(functionals)

    if not iupac.sub_suffix and not functionals:
        res[-1] = lastprime.norm

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
    suffix = all_suffixes2str(iupac)

    return "".join(map(s, [preffix, infix, root, suffix]))


@dataclass
class AtomFeatures:
    chain_index: int
    atom: Atom
    bond_ahead: Bond
    subciupacs: list[IupacName]
    functional_group: GroupMatch | None


@dataclass
class DecFeatures:
    dec: MolDecomposition
    fs: list[AtomFeatures]

    group_scores: list[int] = field(default_factory=list)
    bond_scores: list[int] = field(default_factory=list)
    weak_group_scores: list[int] = field(default_factory=list)
    subchain_scores: list[int] = field(default_factory=list)

    def as_tuple(self):
        return self.dec, self.fs


class Iupac:
    def __init__(self, mol: Molecule):
        self.mol = mol

    @cached_property
    def alg(self):
        return Alogrythms(self.mol)

    @cached_property
    def decomposition(self):
        return self.alg.decompose()

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

    def functional_suffixes(self, features: list[AtomFeatures]):
        result = []
        for f in features:
            if f.functional_group:
                result.append(Synt(f.chain_index, f.functional_group.tag.value))
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
        BondType.SINGLE: 0,
        "weak-functional-group": 1,
        BondType.TRIPLE: 2,
        BondType.DOUBLE: 3,
        "functional-group": 4,
    }

    def feature_bonds(self, feature: AtomFeatures):
        return Iupac.PRIORITIES.get(feature.bond_ahead.type)

    def feature_group(self, feature: AtomFeatures):
        prio = Iupac.PRIORITIES
        if feature.functional_group is None:
            return 0

        group_tag = feature.functional_group.tag.value
        # future:    in weak functional group
        if group_tag in {}:
            return 0

        return prio.get("functional-group")

    def feature_weak_group(self, feature: AtomFeatures):
        prio = Iupac.PRIORITIES
        if feature.functional_group is None:
            return 0

        group_tag = feature.functional_group.tag
        if group_tag not in {}:
            return 0

        return prio.get("weak-functional-group")

    def feature_subchain(self, feature: AtomFeatures):
        alphabetical = 0
        if feature.subciupacs:
            for sub in feature.subciupacs:
                # handle for complex radicals
                alphabetical = max(alphabetical, ord(sub.root_word.norm[0]))

            alphabetical = -alphabetical + ord(Alphabet.MAX.value.norm)
        return alphabetical

    def fpod(self, decomp: MolDecomposition, preffixes: dict[Atom, list[IupacName]]):
        if decomp.is_cycle:
            return self.fpod_cycle(decomp, preffixes)

        straight = decomp
        reversed = decomp.with_chain(tuple(reversechain(straight.chain)))

        return self._fpod_many(preffixes, [straight, reversed])

    def fpod_cycle(
        self, decomp: MolDecomposition, preffixes: dict[Atom, list[IupacName]]
    ):
        features = list(self.features(decomp, preffixes))

        priorities = [
            self.feature_group,
            self.feature_bonds,
            self.feature_weak_group,
            self.feature_subchain,
            # any
            lambda x: 1,
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
        self, preffixes: dict[Atom, list[IupacName]], decs: list[MolDecomposition]
    ):
        features = [
            DecFeatures(dec, list(self.features(dec, preffixes))) for dec in decs
        ]

        #
        for decfs in features:
            decfs.group_scores = list(map(self.feature_group, decfs.fs))

        i = first_max([d.group_scores for d in features])
        if i is not None:
            return features[i].as_tuple()
        #
        for decfs in features:
            decfs.bond_scores = list(map(self.feature_bonds, decfs.fs))

        i = first_max([d.bond_scores for d in features])
        if i is not None:
            return features[i].as_tuple()
        #
        for decfs in features:
            decfs.weak_group_scores = list(map(self.feature_weak_group, decfs.fs))

        i = first_max([d.weak_group_scores for d in features])
        if i is not None:
            return features[i].as_tuple()
        #
        for decfs in features:
            decfs.subchain_scores = list(map(self.feature_subchain, decfs.fs))

        i = first_max([d.subchain_scores for d in features])
        if i is not None:
            return features[i].as_tuple()

        #
        return features[0].as_tuple()

    def decompose_name(self, decomp: MolDecomposition, *, primary: bool = True):
        #
        unordered_preffixes = self.preffixes(decomp)
        #
        decomp, features = self.fpod(decomp, unordered_preffixes)

        preffixes = self.index_preffixes(features, unordered_preffixes)
        root = ROOT_BY_LENGTH[len(decomp.chain)].value
        suffix = self.suffixes_by_features(features, primary=primary)
        func_suffixes = self.functional_suffixes(features)
        infix = Infix.CYCLO.value if decomp.is_cycle else None
        subsuff = Prefix.YL.value if not primary else None

        return IupacName(
            prefixes=preffixes,
            infix=infix,
            root_word=root,
            prime_suffixes=suffix,
            sub_suffix=subsuff,
            func_suffixes=func_suffixes,
            ref=decomp,
        )
