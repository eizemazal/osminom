from dataclasses import dataclass, field
from typing import NamedTuple

from buftinom.algorythms import MolDecomposition
from buftinom.translate import WordForm


def s(obj):
    if not obj:
        return ""
    return str(obj)


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


@dataclass(frozen=True, slots=True, kw_only=True)
class IupacName:
    """
    Structured IUPAC name representation
    """

    subnames: list["IupacName"]
    prefixes: list[Subn] = field(default_factory=list)
    infix: WordForm = None
    func_preffixes: list[Synt]
    root_word: WordForm
    prime_suffixes: list[Synt]
    sub_suffix: Synt = None
    func_suffixes: list[Synt]
    ref: MolDecomposition

    def __repr__(self):
        sn = s(self.subnames)
        if sn:
            sn += " "
        return sn + "".join(
            map(
                s,
                [
                    self.prefixes,
                    self.infix,
                    self.func_preffixes,
                    self.root_word,
                    self.prime_suffixes,
                    self.sub_suffix,
                    self.func_suffixes,
                ],
            )
        )
