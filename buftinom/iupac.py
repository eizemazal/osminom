from collections import defaultdict
from functools import cached_property

from buftinom.algorythms import (
    Alogrythms,
    MolDecomposition,
    SubchainConnection,
)
from buftinom.features import AtomFeatures, Features
from buftinom.lookup import (
    MULTI_BY_PREFIX,
    MULTI_MULTI_BY_PREFIX,
    ROOT_BY_LENGTH,
    Aromatic,
    Infix,
    Prefix,
    PrimarySuffix,
    is_preferred_as_prefix,
    is_preferred_in_prefix,
    provides_split,
)
from buftinom.models import IupacName, Subn, Synt
from buftinom.smileg import Atom, BondType, Molecule
from buftinom.translate import WordForm, WordFormName


def s(obj):
    if not obj:
        return ""
    return str(obj)


def singleidx2str(name: WordForm, idx: int, form: WordFormName):
    if idx > 1:
        return ["-", str(idx), "-", name.get(form)]
    else:
        return [name.get(form)]


def manyidx2str(name: WordForm, ids: list[int], form: WordFormName):
    index = ",".join(map(str, ids))
    multi = MULTI_BY_PREFIX.get(len(ids))
    return ["-", index, "-", multi.value.norm, name.get(form)]


def suffixes2str(suffixes: list[Synt], form: WordFormName):
    """
    Joins the given suffixes to string, selects appropriate word forms,
    Sorts them by chaind indexes, follows IUPAC grammar for indexes and names
    """
    suffixes = [s for s in suffixes if s]
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
    subs, _ = suffixes2str([iupac.sub_suffix], form="norm")
    fform: WordFormName = "norm"
    if subs:
        fform = "sub"

    functionals, _ = suffixes2str(iupac.func_suffixes, form=fform)

    res = primes

    if subs:
        res.extend(subs)
    if functionals:
        res.extend(functionals)

    if lastprime and not iupac.sub_suffix and not functionals:
        res[-1] = lastprime.norm

    return "".join(res)


def func_preffixes2str(iupac: IupacName, *, separate: bool):
    preffixes, _ = suffixes2str(iupac.func_preffixes, "norm")
    preffixes_str = "".join(preffixes)

    if not separate:
        preffixes_str = preffixes_str.removeprefix("-")

    return preffixes_str


def preffixes2str(iupac: IupacName):
    """
    Joins the preffixes of the IupacName to string, selects appropriate word forms,
    Sorts them by chaind indexes, follows IUPAC grammar for indexes and names
    """
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

        name_complex = name[0].isdigit() or " " in name
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


def root2str(iupac: IupacName) -> str:
    if iupac.sub_suffix:
        return iupac.root_word.get("sub")

    if iupac.func_suffixes:
        return iupac.root_word.get("sub")

    return iupac.root_word.get("norm")


def iupac2str(iupac: IupacName) -> str:
    """
    The iupac2str method.

    Converts structural representation of the IupacName to grammatically correct string
    """
    preffix = preffixes2str(iupac)
    infix = iupac.infix
    fpreffix = func_preffixes2str(iupac, separate=bool(preffix or infix))
    root = root2str(iupac)
    suffix = all_suffixes2str(iupac)

    subnames = ""
    if iupac.subnames:
        subnames = " ".join(map(iupac2str, iupac.subnames)) + " "

    return subnames + "".join(map(s, [preffix, infix, fpreffix, root, suffix]))


class Iupac:
    """
    Build Iupac structured representation of the molecule.
    Uses Algorythms class to create decomposition and give names to decomposed parts
    """

    def __init__(self, mol: Molecule):
        self.mol = mol

    @cached_property
    def alg(self):
        return Alogrythms(self.mol)

    @cached_property
    def decomposition(self):
        return self.alg.decompose()

    @cached_property
    def features(self):
        return Features(self.mol)

    def suffixes_by_features(
        self,
        decomp: MolDecomposition,
        features: list[AtomFeatures],
        *,
        primary: bool,
    ):
        """
        Find which primary suffixes should be used based on the bonds of the chain
        """
        result = []

        # aromatic suffix managed by a root word definition
        if decomp.is_aromatic:
            return result

        for feature in features:
            if feature.bond_ahead is None:
                continue

            if feature.bond_ahead.type == BondType.DOUBLE:
                result.append(Synt(feature.chain_index, PrimarySuffix.ENE.value))

            if feature.bond_ahead.type == BondType.TRIPLE:
                result.append(Synt(feature.chain_index, PrimarySuffix.YNE.value))

        if not result and primary:
            result.append(Synt(1, PrimarySuffix.ANE.value))

        return result

    def subsuffix(
        self,
        features: list[AtomFeatures],
        connector: SubchainConnection,
        primary: bool,
    ):
        if primary:
            return None

        for feature in features:
            if not feature.connection:
                continue

            conn = feature.connection

            if conn.via and is_preferred_in_prefix(conn.via.tag):
                return Synt(feature.chain_index, conn.via.tag.value)

            if conn.bond == BondType.SINGLE:
                return Synt(feature.chain_index, Prefix.YL.value)

            if conn.bond == BondType.DOUBLE:
                return Synt(feature.chain_index, Prefix.YLIDENE.value)

            if conn.bond == BondType.TRIPLE:
                return Synt(feature.chain_index, Prefix.YLIDYNE.value)

        raise AssertionError(
            f"Connector {connector} is not connected to the chain {features}"
        )

    def functional_suffixes(self, features: list[AtomFeatures], *, primary: bool):
        """
        Collect all functional suffixed of the chain based on the functional groups it have
        """
        result = []
        for f in features:
            for group in f.functional_groups:
                is_suffix = not is_preferred_in_prefix(group.tag)
                is_suffix = is_suffix and not is_preferred_as_prefix(group.tag)

                if is_suffix:
                    result.append(Synt(f.chain_index, group.tag.value))
        return result

    def children(self, dec: MolDecomposition):
        """
        Recursively find preffixes of the Iupac name.
        """
        subnames: list[IupacName] = []
        preffixes: dict[Atom, list[IupacName]] = defaultdict(list)

        for connection, subchains in dec.connections.items():
            for subchain in subchains:
                #
                if connection.via is None or not provides_split(connection.via.tag):
                    subiupac = self.decompose_name(subchain, primary=False)
                    preffixes[connection.parent].append(subiupac)
                else:
                    subiupac = self.decompose_name(subchain, primary=False)
                    subnames.append(subiupac)
                #

        return subnames, preffixes

    def index_preffixes(
        self, features: list[AtomFeatures], all_preffixes: dict[Atom, list[IupacName]]
    ):
        """
        Give indexes to the found preffixes
        """
        result: list[Subn] = []

        for feature in features:
            if preffixes := all_preffixes.get(feature.atom):
                for preffix in preffixes:
                    result.append(Subn(feature.chain_index, preffix))
        return result

    def functional_preffixes(self, features: list[AtomFeatures]):
        """
        Give indexes to the found functional preffixes
        """
        result: list[Synt] = []

        for feature in features:
            for group in feature.functional_groups:
                if is_preferred_as_prefix(group.tag):
                    result.append(Synt(feature.chain_index, group.tag.value))

        return result

    def root_word(self, decomp: MolDecomposition, features: list[AtomFeatures]):
        if not decomp.is_aromatic:
            return ROOT_BY_LENGTH[len(features)].value

        if len(decomp.chain) == 6:
            return Aromatic.BENZ.value

        raise NotImplementedError(
            f"Aromatic ring with lengh {len(decomp.chain)} is not supported"
        )

    def infix(self, decomp: MolDecomposition, features: list[AtomFeatures]):
        if decomp.is_aromatic or not decomp.is_cycle:
            return None

        return Infix.CYCLO.value

    def decompose_name(self, decomp: MolDecomposition, *, primary: bool = True):
        """
        Recursively convert given decomposition to the IupacName.
        """
        #
        subnames, unordered_preffixes = self.children(decomp)
        #
        decomp, features = self.features.fpod(decomp, unordered_preffixes)

        preffixes = self.index_preffixes(features, unordered_preffixes)
        func_preffixes = self.functional_preffixes(features)
        infix = self.infix(decomp, features)
        #
        root = self.root_word(decomp, features)
        #
        subsuff = self.subsuffix(features, decomp.connected_by, primary)
        suffix = self.suffixes_by_features(decomp, features, primary=primary)
        func_suffixes = self.functional_suffixes(features, primary=primary)

        return IupacName(
            subnames=subnames,
            prefixes=preffixes,
            infix=infix,
            func_preffixes=func_preffixes,
            root_word=root,
            prime_suffixes=suffix,
            sub_suffix=subsuff,
            func_suffixes=func_suffixes,
            ref=decomp,
        )

    def construct_name(self):
        """
        Convert decomposition to the IupacName.
        """
        return self.decompose_name(self.decomposition)
