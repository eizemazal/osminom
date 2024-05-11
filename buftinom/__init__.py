from buftinom.algorythms import Alogrythms
from buftinom.iupac import Iupac, IupacName, iupac2str
from buftinom.smileg import Molecule
from buftinom.smiles_parser import SmilesParser


def _print_debug(iupac: Iupac, iupac_name: IupacName):
    iupac.mol.print_table()
    iupac.decomposition.print()
    print(iupac_name)


def smiles2iupac(smiles: str, *, debug: bool = False) -> list[str]:
    """
    High-level method to translate smiles string to array of IUPAC names associated

    Returns an array of IUPAC names in this smiles string
    """
    mols = SmilesParser(debug=debug).parse(smiles)

    names = []

    for mol in mols:
        iupac = Iupac(mol)
        iupac_name = iupac.construct_name()

        if debug:
            _print_debug(iupac, iupac_name)

        str_name = iupac2str(iupac_name)
        names.append(str_name)

    return names


def smiles2iupac1(smiles: str, *, debug: bool = False) -> str:
    """
    Same as smiles2iupac, but for single molecule

    Returns single string
    """
    names = smiles2iupac(smiles, debug=debug)
    if len(names) != 1:
        raise ValueError(
            f"This method assumes there will be one milecule in smiles, found {names}. "
            + "Consider using smiles2iupac instead."
        )

    return names[0]


__all__ = [
    "Molecule",
    "SmilesParser",
    "Alogrythms",
    "Iupac",
    "IupacName",
    "iupac2str",
    "smiles2iupac",
    "smiles2iupac1",
]
