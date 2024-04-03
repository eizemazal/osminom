from buftinom.lookup import ROOT_BY_LENGTH, PrimarySuffix, bond_prio_cmp
from buftinom.smileg import Atom, BondType, Molecule


class Features:
    def __init__(self, mol: Molecule, chain: list[Atom]):
        self.mol = mol
        self.chain = chain

    def length(self):
        return len(self.chain)

    def fpod(self):
        """first point of difference"""
        for x in zip(
            self.unusual_bonds(reverse=False), self.unusual_bonds(reverse=True)
        ):
            (i, (_, a, _)), (j, (_, b, _)) = x
            if i < j:
                break
            if i > j:
                self.chain = list(reversed(self.chain))
                break
            if i == j:
                bond_prio = bond_prio_cmp(a, b)
                if bond_prio > 0:
                    break
                if bond_prio < 0:
                    self.chain = list(reversed(self.chain))
                    break
                break

    def bonds(self, reverse):
        if len(self.chain) == 1:
            return (0, (self.chain[0], BondType.SINGLE, self.chain[0]))

        chain = self.chain
        if reverse:
            chain = list(reversed(chain))

        for i, (a1, a2) in enumerate(zip(chain, chain[1:])):
            bond = self.mol._bonds.get((a1, a2))
            if not bond:
                raise ValueError(f"Missing bond between {a1} and {a2}")

            yield i, (a1, bond, a2)

    def unusual_bonds(self, reverse=False):
        for i, (a1, bond, a2) in self.bonds(reverse):
            if bond.type != BondType.SINGLE:
                yield i, (a1, bond, a2)


class Iupac:
    def __init__(self, mol: Molecule):
        self.mol = mol

    def suffix_by_bonds(self, bonds: list[tuple[int, tuple[Atom, BondType, Atom]]]):
        if not bonds:
            yield (0, PrimarySuffix.ANE)

        for i, (a1, bond, a2) in bonds:
            match bond.type:
                case BondType.SINGLE:
                    yield (i, PrimarySuffix.ANE)
                case BondType.DOUBLE:
                    yield (i, PrimarySuffix.ENE)
                case BondType.TRIPLE:
                    yield (i, PrimarySuffix.YNE)

    def primary_suffixes(self, primary_suffixes: list[tuple[int, PrimarySuffix]]):
        result = []

        for i, (pos, suffix) in enumerate(primary_suffixes):
            form = suffix.value.short
            if i == len(primary_suffixes) - 1:
                form = suffix.value.norm

            if pos == 0:
                result.append(form)
                continue

            result.extend([str(pos + 1), form])

        return result

    def simple_chain_name(self, chain: list[Atom]):
        f = Features(self.mol, chain)
        f.fpod()

        root = ROOT_BY_LENGTH[f.length()]

        unusual_bonds = list(f.unusual_bonds())
        primary_suffixes = list(self.suffix_by_bonds(unusual_bonds))
        primary_suffixes = sorted(primary_suffixes, key=lambda x: x[1].value.norm)

        psuffixes = self.primary_suffixes(primary_suffixes)
        joined = "-".join(psuffixes)
        link = "-" if psuffixes[0].isnumeric() else ""

        return f"{root.value.norm}{link}{joined}"
