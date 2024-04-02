import ply.lex as lex


class SmilesLexer:
    def __init__(self):
        self.lexer = lex.lex(module=self)

    tokens = (
        "ELEM",
        "LBRACE",
        "RBRACE",
        "BOND",
        "NUM",
        "SLASH",
        "CHIRALITY",
        "CHARGE",
        "SEP",
        "H",
    )

    literals = "()"
    states = (("atom", "inclusive"),)

    def t_SEP(self, t):
        r"\."
        return t

    def t_atom_CHARGE(self, t):
        r"(\+|\-)"
        t.value = 1 if t.value == "+" else -1
        return t

    def t_LBRACE(self, t):
        r"\["
        t.lexer.push_state("atom")
        return t

    def t_atom_RBRACE(self, t):
        r"\]"
        t.lexer.pop_state()
        return t

    def t_atom_CHIRALITY(self, t):
        r"(@@|@)"
        return t

    def t_ANY_SLASH(self, t):
        r"(/|\\)"
        return t

    def t_BOND(self, t):
        r"(-|=|\#|:)"
        return t

    def t_atom_H(self, t):
        r"H"
        return t

    def t_atom_NUM(self, t):
        r"\d+"
        # isotop, nhydrogen, charge, etc.
        t.value = int(t.value)
        return t

    def t_NUM(self, t):
        r"(\d|\%\d+)"
        # Cyclic bond definition
        t.value = t.value
        return t

    def t_ELEM(self, t):
        r"(Br|Cl|B|C|N|O|P|S|F|I|c|n)"
        return t

    def t_atom_ELEM(self, t):
        r"(He|Li|Be|Ne|Na|Mg|Al|Si|Cl|Ar|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu    \
            |Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn \
            |Sb|Te|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu \
            |Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa  \
            |Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|U|P|S|K|Y|I|B|C|N|O|F|c|n)"
        return t

    def t_ANY_error(self, t):
        raise ValueError(f"Invalid token: {t}")

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
