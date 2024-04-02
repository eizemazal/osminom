from __future__ import annotations

import ply.yacc as yacc

from osminom.ast import Ast
from osminom.atom import Atom
from osminom.smiles_lexer import SmilesLexer


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
        """chain : chain MULBOND atom"""
        p[0] = p[1].extend(p[3], p[2])

    def p_chain_brackets(self, p):
        """chain : chain '(' chain ')'"""
        p[0] = p[1].append(p[3])

    def p_chain_brackets_multiple(self, p):
        """chain : chain '(' MULBOND chain ')'"""
        p[0] = p[1].append(p[4], None, p[3])

    def p_chain_brackets_multiple_stereo(self, p):
        """chain : SLASH atom '(' MULBOND atom SLASH chain ')'"""
        if p[4] != "=":
            raise ValueError("Only double bonds are allowed to have stereo")
        bond = "cis=" if p[1] == p[6] else "trans="
        inbra = Ast.from_atom(p[5]).append(p[7])
        p[0] = Ast.from_atom(p[2]).append(inbra, bond=bond)

    # FIXME skip this rule until refactoring conns in Atom
    # def p_chain_mulbond_int(self, p):
    #    """chain : chain MULBOND INT"""
    #    p[0] = [Atom("", connid=p[3]), p[2], p[1]]
    #    return p

    def p_chain_double_stereo(self, p):
        """chain : chain SLASH atom MULBOND atom SLASH"""
        if p[4] != "=":
            raise ValueError("Only double bonds are allowed to have stereo")
        bond = "cis=" if p[2] == p[6] else "trans="
        p[0] = p[1].extend(p[3]).extend(p[5], bond=bond)

    def p_atom(self, p):
        """atom : ELEM
        | atom INT"""
        if len(p) == 2:
            p[0] = Atom(p[1])
        else:
            p[0] = p[1]
            p[0].conns = list(p[2])

    def p_atom_charged(self, p):
        """atom : '[' atom '+' ']'
        | '[' atom  '-' ']'"""
        charge = 1 if p[3] == "+" else -1
        p[0] = p[2]
        p[0].charge = charge

    def p_atom_chiral(self, p):
        "atom : atom CHIRALITY"
        p[0] = p[1]
        p[0].chirality = p[2]

    def p_atom_protonated(self, p):
        """atom : '[' atom 'H' ']'"""
        p[0] = p[2]
        p[0].nprotons = (p[0].nprotons or 0) + 1

    def p_error(self, p):
        raise ValueError("Syntax error")

    def __init__(self, debug=False):
        self.debug = debug
        self.lexer = SmilesLexer()
        self.parser = yacc.yacc(module=self)

    def parse(self, data):
        return self.parser.parse(data, lexer=self.lexer, debug=self.debug)
