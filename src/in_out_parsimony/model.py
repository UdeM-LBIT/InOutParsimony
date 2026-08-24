from ast import literal_eval
from collections.abc import Set
from dataclasses import dataclass

from sowing import Node, traversal
from sowing.repr import newick


@dataclass(frozen=True, slots=True)
class Synteny:
    """Internal divergence node of a synteny tree."""

    # Unique identifier of the node
    name: str

    # Set of genes inside the node
    contents: Set[str] = frozenset()


@dataclass(frozen=True, slots=True)
class Gain(Synteny):
    """Unary node representing the gain of a set of genes."""

    # Set of genes gained from this node to its child
    gained: Set[str] = frozenset()


@dataclass(frozen=True, slots=True)
class Loss(Synteny):
    """Unary node representing the loss of a set of genes."""

    # Set of genes lost from this node to its child
    lost: Set[str] = frozenset()


def read_tree(data: str) -> Node[Synteny, None]:
    """Convert an input tree from Newick to a recursive data structure."""
    return traversal.map(
        lambda node, edge, *_: (
            Synteny(
                name=node.get("name", ""),
                contents=frozenset(literal_eval(node.get("contents", "set()"))),
            ),
            None,
        ),
        traversal.depth(newick.parse(data)),
    )
