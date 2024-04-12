"""
    smiles2iupac script,

    invoke with

    ```bash
    python -m buftinom --help
    ```
"""
import argparse
from dataclasses import dataclass

from buftinom import Iupac, SmilesParser, iupac2str
from buftinom.smileg import debug_atoms


@dataclass
class AppArgs:
    smiles: str
    verbose: bool
    quiet: bool


def parse_args():
    arp = argparse.ArgumentParser(
        "smiles2iupac",
        usage="construct an IUPAC name for given molecule",
    )
    arp.add_argument("smiles", help="molecule in smiles format")
    arp.add_argument(
        "-v",
        "--verbose",
        help="show logs and debug info",
        action="store_true",
    )
    arp.add_argument(
        "-q",
        "--quiet",
        help="print only final result (overrides -v)",
        action="store_true",
    )
    args = AppArgs(**arp.parse_args().__dict__)
    if args.quiet:
        args.verbose = False

    return args


def main():
    args = parse_args()
    debug_atoms(args.verbose)

    parser = SmilesParser(debug=args.verbose)

    molecules = parser.parse(args.smiles)

    if not args.quiet:
        s = "s" if len(molecules) > 1 else ""
        print(f"🔎 Found {len(molecules)} molecule{s}")

    for mol in molecules:
        iupac = Iupac(mol)

        if args.verbose:
            print("🔥 Molecule:")
            mol.print_table()
            print("✂️ Decomposition:")
            iupac.decomposition.print()

        iupac_name = iupac.construct_name()

        if args.verbose:
            print("🏭 Name Construct:")
            print(iupac_name)

        str_name = iupac2str(iupac_name)

        if not args.quiet:
            print(f"🥳 {args.smiles} IUPAC name is:")

        print(str_name)


if __name__ == "__main__":
    main()
