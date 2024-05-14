from __future__ import annotations

import ply.yacc as yacc

from osminom.common.ast import Ast
from osminom.common.atom import Atom
from osminom.common.smiles_lexer import SmilesLexer


class SmilesParser:
    tokens = SmilesLexer.tokens

    def p_smiles(self, p):
        """smiles : chain"""
        p[0] = p[1].pop_root()

    def p_chain_atom(self, p):
        """chain : atom
        | chain atom"""
        if len(p) == 2:
            p[0] = Ast.from_atom(p[1]).push_root()
        elif len(p) == 3:
            p[0] = p[1].extend(p[2])

    def p_chain_atom_multiple_bond(self, p):
        """chain : chain BOND atom"""
        p[0] = p[1].extend(p[3], p[2])

    def p_chain_brackets(self, p):
        """chain : chain '(' chain ')'"""
        p[0] = p[1].append(p[3])

    def p_chain_brackets_multiple(self, p):
        """chain : chain '(' BOND chain ')'"""
        p[0] = p[1].append(p[4], None, p[3])

    def p_chain_brackets_multiple_stereo(self, p):
        """chain : SLASH atom '(' BOND atom SLASH chain ')'"""
        if p[4] != "=":
            raise ValueError("Only double bonds are allowed to have stereo")
        bond = "/" if p[1] == p[6] else "\\"
        inbra = Ast.from_atom(p[5]).append(p[7])
        p[0] = Ast.from_atom(p[2]).append(inbra, bond=bond)

    def p_chain_mulbond_int(self, p):
        """chain : chain BOND INT"""
        # p[0] = [Atom("", closure_id=p[3]), p[2], p[1]]
        p[0] = p[1]
        mapped_id = self.map_closure(p[3])
        p[0].root.atom.closures[mapped_id] = p[2]
        return p

    def p_chain_double_stereo(self, p):
        """chain : chain SLASH atom BOND atom SLASH"""
        if p[4] != "=":
            raise ValueError("Only double bonds are allowed to have stereo")
        bond = "/" if p[2] == p[6] else "\\"
        p[0] = p[1].extend(p[3]).extend(p[5], bond=bond)

    def p_atom(self, p: tuple[None, str]):
        """atom : ELEM
        | atom_bra
        """
        if isinstance(p[1], Atom):
            p[0] = p[1]
        else:
            p[0] = Atom(p[1])

    def p_chain_closure(self, p):
        """chain : chain INT
        | chain '%' INT"""
        p[0] = p[1]
        atom = p[0].root.atom
        if len(p) == 3:
            value = int(p[2])
            while value:
                mapped_id = self.map_closure(value % 10)
                atom.closures[mapped_id] = "-"
                value = value // 10
        else:
            mapped_id = self.map_closure(p[3])
            atom.closures[mapped_id] = "-"

    def p_atom_bra(self, p: tuple[None, str, int, str, int, int, str]):
        """atom_bra : LBRACKET isotope elem_or_h chirality explicit_h atomcharge RBRACKET"""
        _, _, p[3].isotope, _, p[3].chirality, p[3].nprotons, p[3].charge, _ = p
        if p[5] is not None:
            p[3].explicit_protons = True
        p[0] = p[3]

    def p_elem_or_h(self, p):
        """elem_or_h : ELEM
        | H"""
        p[0] = Atom(p[1])

    def p_isotope(self, p: tuple[None, int]):
        """isotope : INT
        |
        """
        if len(p) == 2:
            p[0] = p[1]

    def p_chirality(self, p: tuple[None, str]):
        """chirality : CHIRALITY
        |
        """
        if len(p) == 2:
            p[0] = p[1]

    def p_explicit_h(self, p: tuple[None, int]):
        """explicit_h : H INT
        | H
        |
        """
        p[0] = None
        if len(p) == 2:
            p[0] = 1
        if len(p) == 3:
            p[0] = int(p[2])

    def p_atomcharge(self, p: tuple[None, int]):
        """atomcharge : CHARGE
        | CHARGE INT
        |
        """
        if len(p) == 1:
            p[0] = 0
            return
        p[0] = int(p[1] + str(p[2] if len(p) == 3 else 1))

    def p_atomcharge_charge(self, p: tuple[None, int]):
        """atomcharge : CHARGE CHARGE"""
        # deprecated rule
        if p[1] != p[2]:
            raise ValueError("Invalid charge specification for atom")
        p[0] = int(p[1] + "2")

    def p_error(self, p):
        raise ValueError("Syntax error")

    def map_closure(self, id) -> id:
        if id in self.closure_map:
            return self.closure_map.pop(id)
        else:
            self.closure_map[id] = self.next_closure_id
            self.next_closure_id += 1
            return self.closure_map[id]

    def __init__(self, debug=False):
        self.debug = debug
        self.lexer = SmilesLexer()
        self.parser = yacc.yacc(module=self)
        self.closure_map: dict[int, int] = {}
        self.next_closure_id = 1

    def parse(self, data):
        return self.parser.parse(data, lexer=self.lexer, debug=self.debug)
