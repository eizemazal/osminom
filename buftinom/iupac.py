from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from typing import NamedTuple

from buftinom.algorythms import Alogrythms, FpodReliability, MolDecomposition
from buftinom.lookup import (
    MULTI_BY_PREFIX,
    MULTI_MULTI_BY_PREFIX,
    ROOT_BY_LENGTH,
    Infix,
    Prefix,
    PrimarySuffix,
    bond_prio_cmp,
)
from buftinom.smileg import Atom, Bond, BondType, Molecule
from buftinom.translate import WordForm


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

    def cycle_features(self, decomposition: MolDecomposition):
        chain = decomposition.chain + (decomposition.chain[0],)

        for i, (a1, a2) in enumerate(zip(chain, chain[1:]), start=1):
            bond = self.mol.bonds[(a1, a2)]

            yield AtomFeatures(i, a1, bond, decomposition.connections.get(a1, []))

    def fpod_cycle(self, decomposition: MolDecomposition):
        features = self.cycle_features(decomposition)

        return decomposition, features, FpodReliability.RELIABLE

    def fpod(self, decomposition: MolDecomposition):
        """First point of difference

        Should return correctly oriented decomposition
        """
        # if decomposition.is_cycle:
        #     return self.fpod_cycle(decomposition)

        straight = decomposition
        reverse = MolDecomposition(
            tuple(reversed(straight.chain)),
            straight.connections,
            is_cycle=straight.is_cycle,
        )

        straight_features = list(self.features(straight))
        reverse_features = list(self.features(reverse))

        nsubchains = 0
        nbonds = 0
        rel = FpodReliability.RELIABLE

        subchain_choise = None

        for fstraight, freversed in zip(straight_features, reverse_features):
            # First subchain
            # Do not return immediately, bonds have higher priority
            if subchain_choise is None:
                sub_strait = len(fstraight.subchains)
                sub_revers = len(freversed.subchains)

                subchain_cmp = sub_strait - sub_revers
                if subchain_cmp > 0:
                    subchain_choise = straight, straight_features, rel
                if subchain_cmp < 0:
                    subchain_choise = reverse, reverse_features, rel

                nsubchains += sub_strait

            # First posin of unsaturation
            bond_cmp = bond_prio_cmp(fstraight.bond_ahead, freversed.bond_ahead)
            if bond_cmp > 0:
                return straight, straight_features, rel
            if bond_cmp < 0:
                return reverse, reverse_features, rel

            if fstraight.bond_ahead.type != BondType.SINGLE:
                nbonds += 1

        if subchain_choise is not None:
            return subchain_choise

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
        infix = Infix.CYCLO.value if decomp.is_cycle else None
        subsuff = Prefix.YL.value if not primary else None
        preffix = self.preffixes_by_features(features)

        return IupacName(
            prefixes=preffix,
            infix=infix,
            root_word=root,
            prime_suffixes=suffix,
            sub_suffix=subsuff,
            func_suffix=None,
            ref=decomp,
        )
