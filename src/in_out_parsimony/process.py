from collections.abc import Mapping, Set
from dataclasses import dataclass
from functools import reduce

from sowing import Node, traversal
from sowing.indexed import IndexedTree

from .model import Gain, Loss, Synteny
from .signature import SplitInfo


@dataclass(frozen=True, slots=True)
class LcaMin:
    # Unique identifier of the node
    name: str

    # Minimum set of genes that must be present in the node
    min: Set[str] = frozenset()

    # Minimum set of genes that must be lost relative to the parent node
    min_loss: Set[str] = frozenset()

    # Set of genes for which this node is the lowest common ancestor
    lca: Set[str] = frozenset()

    # Set of genes that have their lowest common ancestor at or below this node
    lca_below: Set[str] = frozenset()


def _union[T](sets: Set[T]) -> Set[T]:
    return reduce(frozenset.__or__, sets, initial=frozenset())


def compute_lca_min(node: Node[Synteny, None]) -> IndexedTree[LcaMin, None]:
    """Preprocess an input tree to compute the minimum and lca sets."""
    contents_below = {}

    for cursor in traversal.depth(node, preorder=False):
        data = cursor.node.data

        if cursor.is_leaf():
            contents_below[data.name] = data.contents
        else:
            contents_below[data.name] = _union(
                contents_below[child.node.data.name] for child in cursor.children()
            )

    lca_below = {}

    for cursor in traversal.depth(node, preorder=True):
        data = cursor.node.data
        lca_below[data.name] = (
            contents_below[data.name]
            if cursor.is_root()
            else lca_below[cursor.up().node.data.name]
            - _union(
                contents_below[sibling.node.data.name] for sibling in cursor.siblings()
            )
        )

    lca = {}
    min_contents = {}
    min_loss = {}
    min_loss[node.data.name] = frozenset()

    for cursor in traversal.depth(node, preorder=False):
        data = cursor.node.data

        lca[data.name] = lca_below[data.name] - _union(
            lca_below[child.node.data.name] for child in cursor.children()
        )

        if cursor.is_leaf():
            min_contents[data.name] = data.contents
        else:
            min_contents_at = _union(
                min_contents[child.node.data.name] - lca[child.node.data.name]
                for child in cursor.children()
            )
            min_contents[data.name] = min_contents_at

            for child in cursor.children():
                child_name = child.node.data.name
                min_loss[child_name] = min_contents_at - min_contents[child_name]

    out_node = traversal.map(
        lambda assoc: LcaMin(
            name=assoc.name,
            lca=lca_below[assoc.name],
            lca_below=lca_below[assoc.name],
            min=min_contents[assoc.name],
            min_loss=min_loss[assoc.name]
        ),
        traversal.depth(node),
    )
    return IndexedTree(out_node)


def _resolve_contents(
    sol, lca_min: IndexedTree[LcaMin, None]
) -> Mapping[str, Set[str]]:
    """
    First postprocessing step to compute the set of contents at each node.

    :param sol: solution tree obtained from the grammar
    :param lca_min: preprocessed input tree obtained from :func:`compute_lca_min`
    :returns: dictionary indicating the contents at each node
    """
    contents = {}

    # Bubble up gains
    for sol_cursor in traversal.depth(sol, preorder=False):
        operator = sol_cursor.node.data.operator
        info = sol_cursor.node.data.args[0]

        if isinstance(info, SplitInfo):
            input_cursor = lca_min[info.name]
            contents[info.name] = input_cursor.node.data.min

            if operator != "leaf":
                left = sol_cursor.down(0).node.data.args[0].name
                right = sol_cursor.down(1).node.data.args[0].name

                in_left = contents[left] - contents[info.name]
                assert len(in_left) >= info.in_left
                contents[info.name] |= frozenset(sorted(in_left)[: info.in_left])

                in_right = contents[right] - contents[info.name]
                assert len(in_right) >= info.in_right
                contents[info.name] |= frozenset(sorted(in_right)[: info.in_right])

    # Push down losses
    for sol_cursor in traversal.depth(sol, preorder=True):
        operator = sol_cursor.node.data.operator
        info = sol_cursor.node.data.args[0]

        if isinstance(info, SplitInfo) and info.out > 0:
            input_cursor = lca_min[info.name]

            assert not input_cursor.is_root()
            up = input_cursor.up().node.data.name

            out = contents[up] - contents[info.name]
            assert len(out) >= info.out
            contents[info.name] |= frozenset(sorted(out)[: info.out])

    return contents


def _resolve_events(
    sol,
    lca_min: IndexedTree[LcaMin, None],
    contents: Mapping[str, Set[str]],
) -> Node[Synteny, None]:
    """
    Second postprocessing step to turn a solution tree into an event tree.

    :param sol: solution tree obtained from the grammar
    :param lca_min: preprocessed input tree obtained from :func:`compute_lca_min`
    :param contents: dictionary obtained from :func:`resolve_contents`
    """

    def present(cursor):
        operator = cursor.node.data.operator
        info = cursor.node.data.args[0]

        match operator:
            case "leaf" | "split":
                event = Synteny(
                    name=info.name,
                    contents=contents.get(info.name),
                )

            case "gain":
                if info.size == 0:
                    return cursor.replace(node=cursor.down().node)

                below_contents = contents[info.name]

                if lca_min[info.name].is_root():
                    above_contents = frozenset()
                else:
                    above_contents = contents[lca_min[info.name].up().node.data.name]

                gained = below_contents - above_contents
                event = Gain(
                    name=info.name,
                    contents=below_contents - gained,
                    gained=gained,
                )

            case "loss":
                if info.size == 0:
                    return cursor.replace(node=cursor.down().node)

                below_contents = contents[info.name]
                above_contents = contents[lca_min[info.name].up().node.data.name]

                lost = above_contents - below_contents
                event = Loss(
                    name=info.name,
                    contents=above_contents,
                    lost=lost,
                )

        return cursor.replace(node=cursor.node.replace(data=event))

    return traversal.fold(present, traversal.depth(sol, preorder=False))


def resolve_solution(sol, lca_min: IndexedTree[LcaMin, None]) -> Node[Synteny, None]:
    """
    Postprocess a solution tree to transform it into an event tree.

    :param sol: solution tree obtained from the grammar
    :param lca_min: preprocessed input tree obtained from :func:`compute_lca_min`
    """
    contents = _resolve_contents(sol, lca_min)
    return _resolve_events(sol, lca_min, contents)


def solution_cost(sol, loss, gain) -> float:
    """
    Compute the cost of a solution under a given cost model.

    :param sol: solution tree obtained from the grammar
    :param loss: loss cost function
    :param gain: gain cost function
    """

    def compute_cost(cursor):
        operator = cursor.node.data.operator
        info = cursor.node.data.args[0]
        value = sum(child.node.data for child in cursor.children())

        match operator:
            case "loss":
                value += loss(info)

            case "gain":
                value += gain(info)

        return cursor.replace(node=Node(data=value))

    return traversal.fold(compute_cost, traversal.depth(sol, preorder=False)).data
