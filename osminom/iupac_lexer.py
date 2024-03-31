import ply.lex as lex
from osminom.ast import Ast
from osminom.smiles_parser import SmilesParser as SP
from copy import deepcopy


aliphatic_roots = ["мет", "эт", "проп", "бут", "пент", "гекс", "гепт", "окт", "нон"]
chains = {
    root: SP().parse("C" * (idx + 1)).label(range(1, idx + 2))
    for idx, root in enumerate(aliphatic_roots)
}


class IupacLexer:

    tokens = (
        "ALIPHATIC",
        "SUFFIX",
        "TBASE",
        "PRFXFORM",
        "TRADICAL",
        # "SUFXFORM",
        # "COMMA",
        "DASH",
        "YL",
        "CYCLO",
        "INT",
    )

    t_ignore = "\t"
    literals = "()[]"

    def t_ALIPHATIC(self, t):
        t.value = deepcopy(chains[t.value])
        return t

    t_ALIPHATIC.__doc__ = "|".join(aliphatic_roots)

    def t_SUFFIX(self, t):
        r"(анол|анамин|ан|ен|ин|аль)"
        values = {
            "анол": SP().parse("O"),
            "ан": Ast(),
            "анамин": SP().parse("N"),
        }
        t.value = values[t.value]
        return t

    def t_TBASE(self, t):
        r"(бензол|карбодиимид)"
        values = {
            "бензол": SP().parse("c1ccccc1").label(range(1, 7)),
            "карбодиимид": SP().parse("N=C=N").label([1, 2, 3]),
        }
        t.value = values[t.value]
        return t

    def t_TRADICAL(self, t):
        r"(аллил|амино|винил|фенил)"
        values = {
            "аллил": SP().parse("CC=C").label(range(1, 4)),
            "амино": SP().parse("N"),
            "винил": SP().parse("C=C").label(range(1, 3)),
            "фенил": SP().parse("c1ccccc1").label(range(1, 7)),
        }
        t.value = values[t.value]
        return t

    def t_PRFXFORM(self, t):
        r"(гидрокси|хлор|циано)"
        values = {
            "гидрокси": SP().parse("O"),
            "хлор": SP().parse("Cl"),
            "циано": SP().parse("C#N"),
        }
        t.value = values[t.value]
        return t

    def t_SUFXFORM(self, t):
        r"(амин|ол|хлорид|цианид)"
        values = {
            "амин": SP().parse("N"),
            "ол": SP().parse("O"),
            "хлорид": SP().parse("Cl"),
            "цианид": SP().parse("C#N"),
        }
        t.value = values[t.value]
        return t

    def t_COMMA(self, t):
        r","
        return t

    def t_DASH(self, t):
        r"-"
        return t

    def t_YL(self, t):
        r"ил"
        return t

    def t_CYCLO(self, t):
        r"цикло"
        return t

    def t_INT(self, t):
        r"\d+"
        t.value = int(t.value)
        return t

    def t_error(self, t):
        print(f"Unexpected symbol {t.value[0]}")
        self.lexer.skip(1)

    def __init__(self, **kwargs):
        self.lexer = lex.lex(module=self, **kwargs)

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
