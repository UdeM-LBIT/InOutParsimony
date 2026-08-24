from collections.abc import Callable
from dataclasses import dataclass

from padapto.signature import Signature


@dataclass(frozen=True, slots=True)
class SplitInfo:
    # Identifier of the corresponding internal node in the input tree
    name: str

    # Number of supplementary genes coming from the left subtree
    in_left: int = 0

    # Number of supplementary genes coming from the right subtree
    in_right: int = 0

    # Number of supplementary genes that are not used in this subtree
    out: int = 0


@dataclass(frozen=True, slots=True)
class BranchInfo:
    # Identifier of the next node on the branch in the input tree
    name: str

    # Number of gained or lost genes on this event
    size: int = 0


@dataclass(frozen=True)
class InOutParsimonySignature[T](Signature[T]):
    """Operators to represent solutions to the in-out parsimony problem."""

    # Leaf of the tree
    leaf: Callable[[SplitInfo], T]

    # Internal node of the tree with two subsolutions
    split: Callable[[SplitInfo, T, T], T]

    # Gain event on a branch of the tree
    gain: Callable[[BranchInfo, T], T]

    # Loss event on a branch of the tree
    loss: Callable[[BranchInfo, T], T]
