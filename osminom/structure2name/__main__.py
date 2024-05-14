"""
    smiles2iupac script,

    invoke with

    ```bash
    python -m buftinom --help
    ```
"""

import argparse
from dataclasses import dataclass

from osminom.common.smiles_parser import SmilesParser
from osminom.structure2name import Iupac, iupac2str
from osminom.structure2name.molecule import debug_atoms
from osminom.structure2name.translate import (
    available_languages,
    get_os_lang,
    override_lang,
)


@dataclass
class AppArgs:
    smiles: str
    verbose: bool
    quiet: bool
    lang: str


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
    arp.add_argument(
        "--lang",
        choices=["os"] + available_languages(),
        default="os",
        help="Translation language, OS language by default",
    )
    args = AppArgs(**arp.parse_args().__dict__)
    if args.quiet:
        args.verbose = False

    if args.lang == "os":
        args.lang = get_os_lang()

    return args


def main(args: AppArgs):
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
    args = parse_args()
    with override_lang(args.lang):
        main(args)
