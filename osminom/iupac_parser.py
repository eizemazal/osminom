import ply.yacc as yacc
from osminom.iupac_lexer import IupacLexer


class IupacParser:

    tokens = IupacLexer.tokens

    def p_name(self, p):
        """name : base"""
        p[0] = p[1]

    def p_base(self, p):
        """base : TBASE
        | chain SUFFIX"""
        if len(p) == 2:
            p[0] = p[1]
        elif len(p) == 3:
            p[0] = p[1].append(p[2], to="1")

    def p_alicyclic(self, p):
        """alicyclic : CYCLO ALIPHATIC"""
        if len(p[2]) < 3:
            raise ValueError("Unable to build cycle of less than 3 atoms")
        p[0] = p[2].cyclize()

    def p_chain(self, p):
        """chain : ALIPHATIC
        | alicyclic"""
        p[0] = p[1]

    def p_base_lform(self, p):
        """base : lform base
        | llocant lform base
        | llocant lform DASH base
        """
        match len(p):
            case 3:
                p[0] = p[2].append(p[1])
            case 4:
                p[0] = p[3].append(p[2], p[1])
            case 5:
                p[0] = p[4].append(p[2], p[1])

    def p_lform(self, p):
        """lform : radical
        | PRFXFORM
        | lform radical
        | llocant lform radical"""
        if len(p) == 2:
            p[0] = p[1]
        if len(p) == 3:
            p[0] = p[2].append(p[1])
        elif len(p) == 4:
            p[0] = p[3].append(p[2], to=p[1])

    def p_lform_brackets(self, p):
        """lform : '(' lform ')'"""
        p[0] = p[2]

    def p_radical_tradical(self, p):
        """radical : TRADICAL"""
        p[0] = p[1]

    def p_radical_chain(self, p):
        """radical : chain YL"""
        p[0] = p[1].rebase(p[1].find(lambda x: x.atom.label == "1"))

    def p_llocant(self, p):
        """llocant : INT DASH"""
        p[0] = int(p[1])

    def p_error(self, p):
        print("Syntax error")

    def __init__(self, debug=False):
        self.debug = debug
        self.lexer = IupacLexer()
        self.parser = yacc.yacc(module=self)

    def parse(self, data):
        return self.parser.parse(data, lexer=self.lexer, debug=self.debug)
