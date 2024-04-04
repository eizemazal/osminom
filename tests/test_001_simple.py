from osminom.iupac_parser import IupacParser
import pytest

# from osminom.iupac_lexer import IupacLexer
# def test_bla(parser):
#    lex = IupacLexer()
#    lex.test("метилметан")
#    exit()


@pytest.fixture
def parser():
    return IupacParser(debug=True)


def test_base_suffix(parser):
    mol = parser.parse("метан")
    assert mol.smiles == "C"
    mol = parser.parse("этан")
    assert mol.smiles == "CC"
    mol = parser.parse("пропан")
    assert mol.smiles == "CCC"
    mol = parser.parse("бутан")
    assert mol.smiles == "CCCC"
    mol = parser.parse("пентан")
    assert mol.smiles == "CCCCC"
    mol = parser.parse("гексан")
    assert mol.smiles == "CCCCCC"
    mol = parser.parse("гептан")
    assert mol.smiles == "CCCCCCC"


def test_base_suffix2(parser):
    mol = parser.parse("октанол")
    assert mol.smiles == "CCCCCCCCO"
    mol = parser.parse("циклопропанамин")
    assert mol.smiles == "C1CC1N"


def test_cyclo(parser):
    mol = parser.parse("циклобутан")
    assert mol.smiles == "C1CCC1"
    with pytest.raises(ValueError):
        mol = parser.parse("циклоэтан")
    with pytest.raises(ValueError):
        mol = parser.parse("циклометан")
    mol = parser.parse("циклогексан")
    assert mol.smiles == "C1CCCCC1"


def test_trivial(parser):
    mol = parser.parse("бензол")
    assert mol.smiles == "c1ccccc1"


def test_prefxform(parser):
    mol = parser.parse("аминобензол")
    assert mol.smiles == "Nc1ccccc1"
    mol = parser.parse("цианоэтан")
    assert mol.smiles == "CCC#N"
    mol = parser.parse("фенилциклопропан")
    assert mol.smiles == "C1CC1c1ccccc1"


def test_radical(parser):
    mol = parser.parse("метилциклогексан")
    assert mol.smiles == "CC1CCCCC1"
    mol = parser.parse("метилметан")
    assert mol.smiles == "CC"
    mol = parser.parse("метилгидроксициклогексан")
    assert mol.smiles == "CC1(O)CCCCC1"
    mol = parser.parse("аллилциклобутан")
    assert mol.smiles == "C1CCC1CC=C"


def test_locant_lform(parser):
    mol = parser.parse("3-гидроксигексанол")
    assert mol.smiles == "CCCC(O)CCO"
    mol = parser.parse("3-(гидрокси)гексанол")
    assert mol.smiles == "CCCC(O)CCO"
    mol = parser.parse("3-(гидрокси)-гексанол")
    assert mol.smiles == "CCCC(O)CCO"
    mol = parser.parse("1-гидрокси-3-метилбензол")
    assert mol.smiles == "c1ccc(C)cc1O"
    mol = parser.parse("(гидроксиэтил)бензол")
    assert mol.smiles == "OC(C)c1ccccc1"
    mol = parser.parse("(2-хлораллил)бензол")
    assert mol.smiles == "ClC(=C)Cc1ccccc1"
    # mol = parser.parse("(3-гидрокси)гексанол")
    # assert mol.smiles == "CCCC(O)CCO"


def test_nested_lform(parser):
    # mol = parser.parse("2-(гидрокси)-2-этилбутанол")
    # assert mol.smiles == "CCC(O)CO"
    mol = parser.parse("2-(2-гидроксиэтил)бутанол")
    assert mol.smiles == "CCC(CO)CCO"

    mol = parser.parse("(гидроксиметил)циклогексан")
    assert mol.smiles == "OCC1CCCCC1"
    mol = parser.parse("(2-(метиламино)этил)бензол")
    assert mol.smiles == "CNCCc1ccccc1"

    # mol = parser.parse("(2-(метиламино)этил)бензол")
    # assert mol.smiles == "CNCCc1ccccc1"

    # mol = parser.parse("(2-(4-гидроксифенил)винил)этанол")
    # assert mol.smiles == "CNCCc1ccccc1"

    # mol = parser.parse("гидроксиэтиламин")
    # assert mol.smiles == "OCCN"


def test_multiprefixes(parser):
    mol = parser.parse("1,2,2-трихлорбутан")
    assert mol.smiles == "CCC(Cl)(Cl)CCl"

    mol = parser.parse("5,5-диэтилциклогексан")
    assert mol.smiles == "CCC(CC)(C1)CCCC1"

    mol = parser.parse("(2,2-(2,4-гидроксифенил)винил)этанол")
    assert mol.smiles == "Oc(c1)cc(O)cc1C(=C)(CCO)c1cc(O)cc(c1)O"

# def test_suffixes(parser):
#    mol = parser.parse("метанол")
#    assert mol.smiles == "CO"
# mol = parser.parse("бутилцианид")
# assert mol.smiles == "CCCCC#N"
