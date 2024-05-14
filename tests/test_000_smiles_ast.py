from osminom.common.ast import Ast, AstNode
from osminom.common.atom import Atom
from osminom.common.smiles_parser import SmilesParser as SP

# import pytest


# @pytest.fixture
# def parser() -> SmilesParser:
#    return SmilesParser(debug=True)


# def test_parser(parser):
#    print(parser.parse("CCCC"))
#    exit()


def test_node_walking_rebase():
    node1 = AstNode(atom=Atom("H"))
    node2 = AstNode(atom=Atom("C"))
    node3 = AstNode(atom=Atom("N"))
    node1.append_child(node2)
    node2.append_child(node3, bond="#")
    assert node2.children == [node3]
    assert node3.parent is node2
    assert node1.children == [node2]
    assert node2.parent is node1
    assert node1.parent is None
    assert not node3.children
    ast = Ast(root=node1)
    assert [node.atom.symbol for node in node1.traverse()] == list("NCH")
    assert [node.upbond for node in node1.traverse()] == list("#--")
    assert [node.atom.symbol for node in node3.traverse_up()] == list("NCH")
    assert [node.upbond for node in node1.traverse()] == list("#--")
    ast.rebase(node3)
    assert node3.children == [node2]
    assert node2.parent is node3
    assert node2.children == [node1]
    assert node1.parent is node2
    assert node3.parent is None
    assert not node1.children
    assert [node.atom.symbol for node in node3.traverse()] == list("HCN")
    assert [node.upbond for node in node3.traverse()] == list("-#-")


def test_smiles():
    mol = SP().parse("CCC")
    assert len(mol) == 3
    assert mol.smiles == "CCC"
    mol = SP().parse("CCCCC")
    assert mol.smiles == "CCCCC"
    mol = SP().parse("c1ccccc1")
    assert len(mol) == 6
    assert mol.smiles == "c1ccccc1"
    mol = SP().parse("CCC#N")
    assert mol.smiles in ("CCC#N", "N#CCC")
    mol = SP().parse("C(C)(C)(C)O")
    assert mol.smiles in ("OC(C)(C)C", "CC(C)(C)O")
    # cumene
    mol = SP().parse("C(C)(C)c1ccccc1O")
    assert mol.smiles == "CC(C)c1ccccc1O"
    # anandamide
    mol = SP().parse("O=C(NCCO)CCCC=CCC=CCC=CCC=CCCCCC")
    assert mol.smiles in (
        "O=C(NCCO)CCCC=CCC=CCC=CCC=CCCCCC",
        "OCCNC(=O)CCCC=CCC=CCC=CCC=CCCCCC",
    )
    # capsaicin
    mol = SP().parse("O=C(NCc1cc(OC)c(O)cc1)CCCCC=CC(C)C")
    assert mol.smiles in (
        "O=C(NCc1cc(OC)c(O)cc1)CCCCC=CC(C)C",
        "c1cc(O)c(OC)cc1CNC(=O)CCCCC=CC(C)C",
    )


def test_traverse():
    mol = SP().parse("CNO")
    assert [node.atom.symbol for node in mol.root.traverse()] == list("ONC")
    mol = SP().parse("C1CCC1")
    assert [node.atom.symbol for node in mol.root.traverse()] == list("CCCC")


def traverse(node: AstNode, indent=0):
    print(" " * indent, node.upbond, node.atom.symbol, sep="")
    for child in node.children:
        traverse(child, indent + 2)


def test_chiral():
    """fix [C@H]"""
    parser = SP(debug=True)
    mol = parser.parse("N[C@H]")

    traverse(mol.root)


def test_mulbond_subchain():
    """fix C(=O)"""
    parser = SP()
    mol = parser.parse("C(=O)")

    traverse(mol.root)
