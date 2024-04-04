from __future__ import annotations

from ast import Tuple

import ply.yacc as yacc

from buftinom.smileg import Atom, Bond, Molecule, MoleculeConstructor
from buftinom.smiles_lexer import SmilesLexer


class SmilesParser:
    def __init__(self, debug=False):
        self.id = 0
        self.debug = debug
        self.lexer = SmilesLexer()
        self.parser = yacc.yacc(module=self)
        self.cosmos: list[MoleculeConstructor] = []

    def next_id(self):
        self.id += 1
        return self.id

    tokens = SmilesLexer.tokens

    def p_smiles(self, p):
        """smiles : molecules"""
        p[0] = p[1]

    def p_molecules_one(self, p):
        """molecules : molecule"""
        p[0] = [p[1]]

    def p_molecules(self, p):
        """molecules : molecules SEP molecule"""
        p[0] = p[1] + [p[3]]

    # Molecule
    #

    def p_molecule(self, p):
        """molecule : atom"""
        _, atom = p
        p[0] = MoleculeConstructor().add_atom(atom)

    def p_molecule_atom(self, p: tuple[MoleculeConstructor, MoleculeConstructor, Atom]):
        """molecule : molecule atom"""
        _, molecule, atom = p
        p[0] = molecule.bind_atom(atom)

    def p_molecule_slash(self, p: tuple[MoleculeConstructor, MoleculeConstructor, str]):
        """molecule : molecule SLASH"""
        _, molecule, slash = p
        p[0] = molecule.push_bond(Bond.make(slash))

    def p_molecule_bond_atom(
        self, p: tuple[MoleculeConstructor, MoleculeConstructor, str]
    ):
        """molecule : molecule BOND"""
        _, molecule, bond = p
        p[0] = molecule.push_bond(Bond.make(bond))

    def p_molecule_closure(
        self, p: tuple[MoleculeConstructor, MoleculeConstructor, str]
    ):
        """molecule : molecule NUM"""
        _, molecule, closure_id = p
        p[0] = molecule.bind_closure(closure_id)

    # Submolecule
    #

    def p_molecule_submolecule(
        self,
        p: tuple[
            MoleculeConstructor, MoleculeConstructor, str, MoleculeConstructor, str
        ],
    ):
        """molecule : molecule '(' submolecule ')'"""
        _, mparent, _, mchild, _ = p
        p[0] = mparent.merge(mchild)

    def p_submolecule(self, p):
        """submolecule : molecule"""
        p[0] = p[1]

    def p_submolecule_bond(
        self, p: tuple[MoleculeConstructor, str, MoleculeConstructor]
    ):
        """submolecule : BOND molecule"""
        _, bond, molecule = p
        molecule.connect_bond = Bond.make(bond)
        p[0] = molecule

    # Atom
    #

    def p_atom(self, p: Tuple[None, str]):
        """atom : atom_open
        | atom_parametrized
        """
        p[0] = p[1]

    def p_atom_open(self, p: tuple[None, str]):
        """atom_open : ELEM"""
        p[0] = Atom(id=self.next_id(), symbol=p[1])

    def p_atom_parametrized(self, p: tuple[None, str, int, str, int, int, str]):
        """atom_parametrized : LBRACE isotope ELEM chirality hydrogen charge RBRACE"""
        _, _, isotope, elem, chirality, hydrogen, charge, _ = p
        p[0] = Atom(
            id=self.next_id(),
            symbol=elem,
            isotope=isotope,
            chirality=chirality,
            hydrogen=hydrogen,
            charge=charge,
        )

    def p_atom_parametrized_h(self, p: tuple[None, str, int, str, int, int, str]):
        """atom_parametrized : LBRACE isotope H chirality charge RBRACE"""
        _, _, isotope, elem, chirality, charge, _ = p
        p[0] = Atom(
            id=self.next_id(),
            symbol=elem,
            isotope=isotope,
            chirality=chirality,
            charge=charge,
        )

    def p_isotope(self, p: tuple[None, int]):
        """isotope : NUM
        | empty
        """
        if len(p) == 2:
            p[0] = p[1]

    def p_chirality(self, p: tuple[None, str]):
        """chirality : CHIRALITY
        | empty
        """
        if len(p) == 2:
            p[0] = p[1]

    def p_hydrogen(self, p: tuple[None, int]):
        """hydrogen : H NUM
        | H
        """
        if len(p) == 2:
            p[0] = 1
        if len(p) == 3:
            p[0] = int(p[2])

    def p_hydrogen_empty(self, p: tuple[None]):
        """hydrogen : empty"""
        pass

    def p_charge(self, p: tuple[None, int]):
        """charge : charge_num
        | charge_chain
        | empty
        """
        if len(p) == 2:
            p[0] = p[1]

    def p_charge_num(self, p: tuple[None, int, int]):
        """charge_num : CHARGE NUM"""
        p[0] = p[1] * p[2]

    def p_charge_chain(self, p: tuple[None, int]):
        """charge_chain : charge_chain CHARGE
        | CHARGE
        """
        if len(p) == 2:
            p[0] = p[1]
        if len(p) == 3:
            p[0] = p[1] + p[2]

    ###

    def p_empty(self, _):
        """empty :"""
        pass

    def p_error(self, p):
        raise ValueError(f"Syntax error: {p}")

    def parse(self, data) -> list[Molecule]:
        return self.parser.parse(data, lexer=self.lexer, debug=self.debug)  #
