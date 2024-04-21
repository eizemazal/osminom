from enum import Enum

from buftinom.smileg import Bond, BondType
from buftinom.translate import translate as _


class RootWord(Enum):
    METH = _("root.meth")
    ETH = _("root.eth")
    PROP = _("root.prop")
    BUT = _("root.but")
    PENT = _("root.pent")
    HEX = _("root.hex")
    HEPT = _("root.hept")
    OCT = _("root.oct")
    NON = _("root.non")
    DEC = _("root.dec")
    UNDEC = _("root.undec")
    DODEC = _("root.dodec")
    TRIDEC = _("root.tridec")
    TETRADEC = _("root.tetradec")
    PENTADEC = _("root.pentadec")
    HEXADEC = _("root.hexadec")
    HEPTADEC = _("root.heptadec")
    OCTADEC = _("root.octadec")
    NONADEC = _("root.nonadec")
    ICOS = _("root.icos")


class PrimarySuffix(Enum):
    # all C-C bonds
    ANE = _("psuffix.ane", short=True)
    # one C=C
    ENE = _("psuffix.ene", short=True)
    # two C=C
    DIENE = _("psuffix.diene", short=True)
    # one C≡C
    YNE = _("psuffix.yne", short=True)
    # two C≡C
    DIYNE = _("psuffix.diyne", short=True)
    # one C=C & one C≡C
    ENYNE = _("psuffix.enyne", short=True)


class FunctionalGroup(Enum):
    CARBOXYLIC_ACID = _("functional-group.carboxylic-acid", sub=True)
    AMIDE = _("functional-group.amide", sub=True)
    NITRILE = _("functional-group.nitrile", sub=True)
    ALDEHYDE = _("functional-group.al")
    KETONE = _("functional-group.one")
    OXY = _("functional-group.oxy")
    AMINO = _("functional-group.amino")
    ALCOHOL = _("functional-group.ol")
    THIOL = _("functional-group.thiol")
    AMINE = _("functional-group.amine")
    IMINE = _("functional-group.imine")
    ESTER = _("functional-group.ester", sub=True)
    BROM = _("functional-group.brom")
    CHLOR = _("functional-group.chlor")


class Infix(Enum):
    CYCLO = _("infix.cyclo")
    SPIRO = _("infix.spiro")
    BICYCLO = _("infix.bicyclo")


class Aromatic(Enum):
    BENZ = _("aromatic.benzene", sub=True)


class Prefix(Enum):
    YL = _("prefix.yl")
    YLIDENE = _("prefix.ylidene")
    YLIDYNE = _("prefix.ylidyne")


class MultiplyingPrefix(Enum):
    ONE = _("empty")
    DI = _("multiplying.di")
    TRI = _("multiplying.tri")
    TETRA = _("multiplying.tetra")
    PENTA = _("multiplying.penta")

    BIS = _("multiplying.bis")
    TRIS = _("multiplying.tris")
    TETRAKIS = _("multiplying.tetrakis")
    PENTAKIS = _("multiplying.pentakis")


class Alphabet(Enum):
    MIN = _("alphabet.min")
    MAX = _("alphabet.max")


BOND_PRIORITY = {
    BondType.SINGLE: 10,
    BondType.DOUBLE: 30,
    BondType.TRIPLE: 20,
}


def bond_prio_cmp(a: Bond, b: Bond) -> int:
    return BOND_PRIORITY[a.type] - BOND_PRIORITY[b.type]


ROOT_BY_LENGTH = {
    1: RootWord.METH,
    2: RootWord.ETH,
    3: RootWord.PROP,
    4: RootWord.BUT,
    5: RootWord.PENT,
    6: RootWord.HEX,
    7: RootWord.HEPT,
    8: RootWord.OCT,
    9: RootWord.NON,
    10: RootWord.DEC,
    11: RootWord.UNDEC,
    12: RootWord.DODEC,
    13: RootWord.TRIDEC,
    14: RootWord.TETRADEC,
    15: RootWord.PENTADEC,
    16: RootWord.HEXADEC,
    17: RootWord.HEPTADEC,
    18: RootWord.OCTADEC,
    19: RootWord.NONADEC,
    20: RootWord.ICOS,
}


MULTI_BY_PREFIX = {
    1: MultiplyingPrefix.ONE,
    2: MultiplyingPrefix.DI,
    3: MultiplyingPrefix.TRI,
    4: MultiplyingPrefix.TETRA,
    5: MultiplyingPrefix.PENTA,
}


MULTI_MULTI_BY_PREFIX = {
    1: MultiplyingPrefix.ONE,
    2: MultiplyingPrefix.BIS,
    3: MultiplyingPrefix.TRIS,
    4: MultiplyingPrefix.TETRAKIS,
    5: MultiplyingPrefix.PENTAKIS,
}


PREFFERED_AS_PREFFIXES = {
    FunctionalGroup.BROM,
    FunctionalGroup.CHLOR,
}


PREFFERED_IN_PREFFIXES = {
    FunctionalGroup.OXY,
    FunctionalGroup.AMINO,
}

CAN_SPLIT_NAME = {
    FunctionalGroup.OXY,
    FunctionalGroup.AMINO,
}


def is_preferred_in_prefix(name: FunctionalGroup):
    if name in PREFFERED_IN_PREFFIXES:
        return True

    return False


def is_preferred_as_prefix(name: FunctionalGroup):
    if name in PREFFERED_AS_PREFFIXES:
        return True

    return False


def provides_split(name: FunctionalGroup):
    if name in CAN_SPLIT_NAME:
        return False

    return True
