from __future__ import annotations
from typing import Iterator
from dataclasses import dataclass, field
from osminom.atom import Atom
from copy import deepcopy

@dataclass
class AstNode:
    atom: Atom
    upbond: str = "-"
    children: list[AstNode] = field(default_factory=list)
    _parent: AstNode | None = None

    @property
    def parent(self) -> AstNode:
        return self._parent

    def get_root(self) -> AstNode:
        r = self
        next2go = self.parent
        while next2go and r in next2go.children:
            r = next2go
            next2go = r.parent
        return r

    def append_child(self, child: AstNode, bond: str = "-") -> AstNode:
        if child._parent is not None:
            raise ValueError("AST inconsistent")
        self.children.append(child)
        child._parent = self
        child.upbond = bond
        return self

    def traverse(self) -> Iterator[AstNode]:
        for child in reversed(self.children):
            for r in child.traverse():
                yield r
        yield self

    def traverse_up(self) -> Iterator[AstNode]:
        yield self
        if self._parent:
            for r in self._parent.traverse_up():
                yield r

    def to_smiles(self) -> str:
        r = ""
        for idx, child in enumerate(self.children):
            child_smiles = child.to_smiles()
            if child.upbond != "-":
                child_smiles = child.upbond + child_smiles
            if idx < len(self.children) - 1:
                child_smiles = "(" + child_smiles + ")"
            r += child_smiles
        return self.atom.smiles_repr() + r


@dataclass
class Ast:
    root: AstNode | None = None
    root_stack: list[AstNode] = field(default_factory=list)
    last_connid: int = 0

    @staticmethod
    def from_atom(atom: Atom) -> Ast:
        return Ast(root=AstNode(atom=atom))

    def __len__(self) -> int:
        return sum([1 for _ in self.root.traverse()])

    def append(self, rhs: Ast, to: str | None = None, bond: str = "-") -> Ast:
        if rhs.root is None:
            return self
        node = self.root if to is None else self.find(lambda x: x.atom.label == str(to))
        if not node:
            raise ValueError(f"Label {to} undefined for the parent structure.")
        if rhs.root_stack:
            rhs.rebase(rhs.root_stack[0])
            rhs.root = rhs.root_stack[0]
        rhs.clearlabels()
        node.append_child(rhs.root, bond)
        return self

    def append_multiple(self, rhs_base: Ast, to: list[str] | None = None, bond: str = "-") -> Ast:
        if rhs_base.root is None:
            return self
        nodes = [self.root] if to is None else self.find_multiple_labels(to)
        # print(nodes)
        if not nodes:
            raise ValueError(f"Label {to} undefined for the parent structure.")

        for node in nodes:
            rhs = deepcopy(rhs_base)
            if rhs.root_stack:
                rhs.rebase(rhs.root_stack[0])
                rhs.root = rhs.root_stack[0]
            rhs.clearlabels()
            node.append_child(rhs.root, bond)
        return self

    def extend(self, atom: Atom, bond: str = "-") -> Ast:
        node = AstNode(atom=atom)
        node.append_child(self.root)
        self.root.upbond = bond
        self.root = node
        return self

    def push_root(self, node: AstNode = None):
        self.root_stack.append(node if node is not None else self.root)
        return self

    def pop_root(self):
        self.root = self.root_stack.pop()
        self.rebase(self.root)
        return self

    def find(self, selector_fn: callable[[AstNode], bool]) -> AstNode | None:
        for node in self.root.traverse():
            if selector_fn(node):
                return node
        return None

    def find_all(self, selector_fn: callable[[AstNode], bool]) -> list[AstNode]:
        return [node for node in self.root.traverse() if selector_fn(node)]

    def find_multiple_labels(self, labels: list[str]) -> list[AstNode]:
        return [self.find(lambda x: x.atom.label == str(label)) for label in labels]

    def rebase(self, node: AstNode) -> Ast:
        nodes_to_root = list(node.traverse_up())
        last = node
        last_upbond = node.upbond
        node.upbond = "-"
        node._parent = None
        for predecessor in nodes_to_root[1:]:
            predecessor._parent = last
            predecessor.children.remove(last)
            last.children.append(predecessor)
            predecessor.upbond, last_upbond = last_upbond, predecessor.upbond
            # last.upbond, predecessor.upbond = predecessor.upbond, last.upbond
            last = predecessor
        self.root = node
        return self

    def label(self, numbers: list[int | str]) -> Ast:
        for node, number in zip(self.root.traverse(), numbers):
            node.atom.label = str(number)
        return self

    def clearlabels(self) -> Ast:
        for node in self.root.traverse():
            node.atom.label = ""

    def cyclize(self) -> Ast:
        self.last_connid += 1
        one = self.root
        another = self.find(lambda x: x.atom.label == "1")
        one.atom.add_conn(self.last_connid)
        another.atom.add_conn(self.last_connid)
        return self

    def _rebase_to_farthest(self) -> None:
        root_distances: tuple[AstNode, int] = [
            (node, sum([1 for _ in node.traverse_up()]))
            for node in self.find_all(lambda x: not x.children)
        ]
        root_distances = sorted(root_distances, key=lambda x: x[1], reverse=True)
        if root_distances:
            self.rebase(root_distances[0][0])

    def canonicalize(self) -> Ast:
        self._rebase_to_farthest()
        self._rebase_to_farthest()
        return self

    @property
    def smiles(self):
        self.canonicalize()
        return self.root.to_smiles()
