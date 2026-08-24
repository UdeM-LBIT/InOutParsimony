
from padapto.structure.grammar import clause, grammar, predicate
from padapto.structure.pattern import (
    Empty,
    Item,
    Term,
    Tree,
    Var,
    Zero,
    chain,
)
from sowing import Node

from .signature import BranchInfo, InOutParsimonySignature, SplitInfo


@grammar
class InOutParsimonyGrammar:
    """Dynamic programming grammar to solve the in-out-parsimony problem."""

    alg: InOutParsimonySignature

    @predicate
    @staticmethod
    def node(tree: Node, gains: int, out: int):
        return  # type: ignore

    @predicate
    @staticmethod
    def branch(tree: Node, gains: int, out: int):
        return  # type: ignore

    @predicate
    @staticmethod
    def branch_gain(tree: Node, gains: int, out: int):
        return  # type: ignore

    @clause(
        predicate="node",
        tree=Tree(
            node=Var("data"),
            children=chain(
                Item(Var("left")),
                Item(Var("right")),
            ),
        ),
        gains=chain(Term(Var("gains_left")), Term(Var("gains_right"))),
        out=Var("out"),
    )
    def _split(self, data, left, right, gains_left, gains_right, out):
        in_left = len(left.node.data.lca_below) - gains_left
        in_right = len(right.node.data.lca_below) - gains_right

        if in_left >= 0 and in_right >= 0:
            return self.alg.split(
                SplitInfo(
                    name=data.name,
                    in_left=in_left,
                    in_right=in_right,
                    out=out,
                ),
                self.branch(
                    tree=left,
                    gains=gains_left,
                    out=out + in_right + len(data.min - left.node.data.min),
                ),
                self.branch(
                    tree=right,
                    gains=gains_right,
                    out=out + in_left + len(data.min - right.node.data.min),
                ),
            )
        else:
            return self.alg.null()

    @clause(
        predicate="node",
        tree=Tree(node=Var("data"), children=Empty()),
        gains=Zero(),
        out=Zero(),
    )
    def _leaf(self, data):
        return self.alg.leaf(SplitInfo(name=data.name))

    @clause(
        predicate="branch",
        tree=Var("tree"),
        gains=Var("gains"),
        out=chain(Term(Var("losses")), Term(Var("out"))),
    )
    def _loss(self, tree, gains, losses, out):
        return self.alg.loss(
            BranchInfo(name=tree.node.data.name, size=losses),
            self.branch_gain(tree=tree, gains=gains, out=out),
        )

    @clause(
        predicate="branch_gain",
        tree=Var("tree"),
        gains=chain(Term(Var("gains")), Term(Var("gains_below"))),
        out=Var("out"),
    )
    def _gain(self, tree, gains, gains_below, out):
        return self.alg.gain(
            BranchInfo(name=tree.node.data.name, size=gains),
            self.node(tree=tree, gains=gains_below, out=out),
        )
