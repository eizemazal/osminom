import ply.lex as lex


class SmilesLexer:
    tokens = (
        "ELEM",
        "INT",
        "MULBOND",
        "SLASH",
        "CHIRALITY",
    )
    literals = "()[]+-H."

    def t_CHIRALITY(self, t):
        r"(@@|@)"
        return t

    def t_SLASH(self, t):
        r"(/|\\)"
        return t

    def t_MULBOND(self, t):
        r"(=|\#)"
        return t

    def t_INT(self, t):
        r"\d+"
        t.value = t.value
        return t

    def t_ELEM(self, t):
        r"(He|Li|Be|Ne|Na|Mg|Al|Si|Cl|Ar|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu    \
            |Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn \
            |Sb|Te|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu \
            |Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa  \
            |Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|U|P|S|K|Y|I|B|C|N|O|F|c|n)"
        return t

    def t_error(self, _):
        raise ValueError("Invalid token")

    def __init__(self):
        self.lexer = lex.lex(module=self)

    def input(self, data):
        return self.lexer.input(data)

    def token(self):
        return self.lexer.token()

    def test(self, data):
        self.lexer.input(data)
        while True:
            tok = self.lexer.token()
            if not tok:
                break
            print(tok)
