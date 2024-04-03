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
    # (all C-C bonds)
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


BOND_PRIORITY = {
    BondType.SINGLE: 1,
    BondType.DOUBLE: 3,
    BondType.TRIPLE: 2,
}


def bond_prio_cmp(a: Bond, b: Bond) -> int:
    return BOND_PRIORITY[a.type] - BOND_PRIORITY[b.type]


class FunctionalGroup(Enum):
    pass


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
