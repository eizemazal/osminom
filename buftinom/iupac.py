from collections import defaultdict, deque
from functools import cached_property

from buftinom.algorythms import (
    Alogrythms,
    MolDecomposition,
    SubchainConnection,
)
from buftinom.features import AtomFeatures, Features
from buftinom.funcgroups import GroupMatch
from buftinom.lookup import (
    MULTI_BY_PREFIX,
    MULTI_MULTI_BY_PREFIX,
    ROOT_BY_LENGTH,
    UNIQUE_SUBSUFFIX,
    Aromatic,
    Infix,
    Prefix,
    PrimarySuffix,
    is_always_preferred_as_prefix,
    is_preferred_as_prefix_in_subchain,
    is_preferred_in_subprefix,
    is_preffix_in_arimatic,
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


def select_form(
    iupac: IupacName,
    groups: list[tuple[WordForm, list[int]]],
    form: WordFormName,
):
    word_end = not iupac.sub_suffixes and not iupac.func_suffixes
    complex_suff = iupac.func_suffixes and (
        len(iupac.func_suffixes) > 1 or iupac.func_suffixes[0].index > 1
    )

    for id, (name, ids) in enumerate(groups):
        if form != "short":
            yield form, name, ids
            continue

        if word_end or complex_suff:
            yield "norm", name, ids
            continue

        if id + 1 >= len(groups):
            yield form, name, ids
            continue

        _, next_ids = groups[id + 1]

        if len(next_ids) > 1 or next_ids[0] > 1:
            yield "norm", name, ids
            continue

        yield form, name, ids


def suffixes2str(iupac: IupacName, suffixes: list[Synt], form: WordFormName, sort=True):
    """
    Joins the given suffixes to string, selects appropriate word forms,
    Sorts them by chaind indexes, follows IUPAC grammar for indexes and names
    """
    suffixes = suffixes and [s for s in suffixes if s]
    if not suffixes:
        return []

    res: list[str] = []

    groupped: dict[WordForm, list[int]] = defaultdict(list)
    for idx, name in suffixes:
        groupped[name].append(idx)

    groups = groupped.items()
    if sort:
        groups = sorted(groupped.items(), key=lambda g: g[0])

    for frm, name, ids in select_form(iupac, groups, form):
        if len(ids) > 1:
            res.extend(manyidx2str(name, ids, frm))
        else:
            res.extend(singleidx2str(name, ids[0], frm))

    return res


def all_suffixes2str(iupac: IupacName):
    primes = suffixes2str(iupac, iupac.prime_suffixes, form="short")
    subs = suffixes2str(iupac, iupac.sub_suffixes, form="norm", sort=False)
    fform: WordFormName = "norm"
    if subs:
        fform = "sub"

    functionals = suffixes2str(iupac, iupac.func_suffixes, form=fform)

    res = primes

    if subs:
        res.extend(subs)
    if functionals:
        res.extend(functionals)

    return "".join(res)


def func_preffixes2str(iupac: IupacName, *, separate: bool):
    preffixes = suffixes2str(iupac, iupac.func_preffixes, "pref")
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
    if iupac.sub_suffixes:
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

    def subsuffixes(
        self,
        features: list[AtomFeatures],
        connector: SubchainConnection,
        *,
        primary: bool,
    ):
        if primary:
            return None

        res = deque([])

        for feature in features:
            if not feature.connection:
                continue

            conn = feature.connection

            if conn.via and is_preferred_in_subprefix(conn.via.tag):
                synt = Synt(feature.chain_index, conn.via.tag.value)
                if conn.via.tag in UNIQUE_SUBSUFFIX:
                    return [synt]

                res.append(synt)

            if conn.bond == BondType.SINGLE:
                res.appendleft(Synt(feature.chain_index, Prefix.YL.value))

            if conn.bond == BondType.DOUBLE:
                res.appendleft(Synt(feature.chain_index, Prefix.YLIDENE.value))

            if conn.bond == BondType.TRIPLE:
                res.appendleft(Synt(feature.chain_index, Prefix.YLIDYNE.value))

        if res:
            return list(res)

        raise AssertionError(
            f"Connector {connector} is not connected to the chain {features}"
        )

    def preffered_in_preffix(
        self, group: GroupMatch, *, primary: bool, decomp: MolDecomposition
    ):
        always = is_always_preferred_as_prefix(group.tag)
        in_subchain = not primary and is_preferred_as_prefix_in_subchain(group.tag)
        aromatic = decomp.is_aromatic and is_preffix_in_arimatic(group.tag)
        flat_flexed = not decomp.is_cycle and group.root_flexed

        return always or in_subchain or aromatic or flat_flexed

    def functional_preffixes(
        self,
        decomp: MolDecomposition,
        features: list[AtomFeatures],
        *,
        primary: bool,
    ):
        """
        Give indexes to the found functional preffixes
        """
        result: list[Synt] = []

        for feature in features:
            for group in feature.functional_groups:
                if is_preferred_in_subprefix(group.tag):
                    continue

                if self.preffered_in_preffix(group, primary=primary, decomp=decomp):
                    result.append(Synt(feature.chain_index, group.tag.value))

        return result

    def functional_suffixes(
        self,
        decomp: MolDecomposition,
        features: list[AtomFeatures],
        *,
        primary: bool,
    ):
        """
        Collect all functional suffixed of the chain based on the functional groups it have
        """
        main = []
        for f in features:
            for group in f.functional_groups:
                if is_preferred_in_subprefix(group.tag):
                    continue

                if not self.preffered_in_preffix(group, primary=primary, decomp=decomp):
                    main.append(Synt(f.chain_index, group.tag.value))

        return main

    def children(self, dec: MolDecomposition):
        """
        Recursively find preffixes of the Iupac name.
        """
        subnames: list[IupacName] = []
        preffixes: dict[Atom, list[IupacName]] = defaultdict(list)

        for connection, subchains in dec.connections.items():
            for subchain in subchains:
                #
                subiupac = self.decompose_name(subchain, primary=False)

                if connection.via is None or not provides_split(connection.via.tag):
                    preffixes[connection.parent].append(subiupac)
                else:
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
        subnames, subname_preffixes = self.children(decomp)
        #
        decomp, features = self.features.fpod(decomp, subname_preffixes)

        preffixes = self.index_preffixes(features, subname_preffixes)
        #
        func_preffixes = self.functional_preffixes(decomp, features, primary=primary)
        infix = self.infix(decomp, features)
        #
        root = self.root_word(decomp, features)
        #
        subsuffixes = self.subsuffixes(features, decomp.connected_by, primary=primary)
        suffix = self.suffixes_by_features(decomp, features, primary=primary)
        func_suffixes = self.functional_suffixes(decomp, features, primary=primary)

        return IupacName(
            subnames=subnames,
            prefixes=preffixes,
            infix=infix,
            func_preffixes=func_preffixes,
            root_word=root,
            prime_suffixes=suffix,
            sub_suffixes=subsuffixes,
            func_suffixes=func_suffixes,
            ref=decomp,
        )

    def construct_name(self):
        """
        Convert decomposition to the IupacName.
        """
        return self.decompose_name(self.decomposition)
