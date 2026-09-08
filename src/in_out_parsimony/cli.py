import itertools
import sys
from argparse import ArgumentParser
from collections import defaultdict
from math import exp
from random import Random

from padapto.circuit import enumerate_solutions, eval_inside, eval_outside, sample
from padapto.evaluation import (
    additive,
    boltzmann,
    count,
    join,
    lex,
    pareto,
    power,
    trace,
)
from padapto.math import logaddexp
from sowing import traversal

from .grammar import InOutParsimonyGrammar
from .model import read_tree
from .process import compute_lca_min, resolve_solution, solution_cost
from .signature import InOutParsimonySignature


def run_algorithm(algebra, lca_min):
    """Run the InOutParsimony algorithm on the given input with the given algebra."""
    gram = InOutParsimonyGrammar(algebra)
    return gram.branch(
        tree=lca_min.root.unzip(),
        gains=len(lca_min.root.data.lca_below),
        out=0,
    )


def do_pareto(lca_min):
    """Compute the set of Pareto-optimal vectors for an InOutParsimony instance."""
    counter = count(InOutParsimonySignature)
    min_gain_open = additive(
        InOutParsimonySignature, gain=lambda event: 1 if event.size > 0 else 0
    )
    min_loss_open = additive(
        InOutParsimonySignature, loss=lambda event: 1 if event.size > 0 else 0
    )
    min_loss_ext = additive(
        InOutParsimonySignature,
        loss=lambda event: event.size,
    )
    count_vecs = (
        join(
            gain_open=min_gain_open,
            loss_open=min_loss_open,
            loss_ext=min_loss_ext,
            count=counter,
        )
        | power()
        | pareto("gain_open", "loss_open", "loss_ext")
    )

    print("# (gain_open, loss_open, loss_ext): count")

    for rec in run_algorithm(count_vecs, lca_min):
        print(
            "("
            + ", ".join(map(str, (rec.gain_open, rec.loss_open, rec.loss_ext)))
            + ")",
            rec.count,
        )


def do_solutions_show(
    lca_min,
    costs,
    temperature,
    gen,
    sample_size,
):
    """Extract and present solutions for an InOutParsimony instance."""
    counter = count(InOutParsimonySignature)
    tracer = trace(InOutParsimonySignature)

    if sample_size == float("inf"):
        sample_range = itertools.count()
    else:
        sample_range = range(int(sample_size))

    if temperature in (0, float("inf")) or sample_size in (0, float("inf")):
        if temperature == 0:
            min_cost = additive(InOutParsimonySignature, **costs)
            tracer_min = join(cost=min_cost, count=counter, circ=tracer) | lex("cost")
            result = run_algorithm(tracer_min, lca_min)
        else:
            tracer_all = join(count=counter, circ=tracer)
            result = run_algorithm(tracer_all, lca_min)

        print(f"# count={result.count}")

        for i, sol in zip(sample_range, enumerate_solutions(result.circ), strict=False):
            print(f"# solution no. {i + 1}, cost={solution_cost(sol, **costs)}")
            print(resolve_solution(sol, lca_min))
    else:
        circ = run_algorithm(tracer, lca_min)
        set_boltzmann = boltzmann(
            InOutParsimonySignature, temperature=temperature, **costs
        )

        for i in sample_range:
            sol = sample(circ, gen, set_boltzmann, log_weights=True)
            print(f"# solution no. {i + 1}, cost={solution_cost(sol, **costs)}")
            print(resolve_solution(sol, lca_min))


def do_solutions_tables(lca_min, costs):
    """Print the dynamic programming tables for an InOutParsimony instance."""
    min_cost = additive(InOutParsimonySignature, **costs)
    gram = InOutParsimonyGrammar(min_cost)

    gene_count = len(lca_min.root.data.lca_below)
    result = gram.branch(tree=lca_min.root.unzip(), gains=gene_count, out=0)

    tables = defaultdict(dict)

    for key, value in gram.memo["node"].items():
        node = key["tree"].node.data.name
        gains = key["gains"]
        out = key["out"]
        tables[node][(gains, out)] = value

    width = max(len(str(result)), 3)

    for node in lca_min:
        print(f"# node={node}")
        print(f"{r'o\g':>{width}}", end=" ")

        for gains in range(gene_count + 1):
            print(f"{gains:>{width}}", end=" ")

        print()

        for out in range(gene_count + 1):
            print(f"{out:>{width}}", end=" ")
            for gains in range(gene_count + 1):
                value = tables.get(node, {}).get((gains, out), float("inf"))
                print(f"{value:>{width}}", end=" ")

            print()

        print()


def do_solutions_avg(lca_min, costs, temperature):
    """Compute the average content size at each node for an InOutParsimony instance."""
    if temperature == 0:
        temperature = sys.float_info.epsilon

    tracer = trace(InOutParsimonySignature)
    circ = run_algorithm(tracer, lca_min)
    set_boltzmann = boltzmann(InOutParsimonySignature, temperature=temperature, **costs)

    boltz_in = eval_inside(circ, set_boltzmann)
    boltz_out = eval_outside(circ, set_boltzmann, boltz_in, log_weights=True)

    class ExpectDict:
        def __init__(self):
            self.values = {}

        def add(self, key, value):
            if key not in self.values:
                self.values[key] = float("-inf")

            self.values[key] = logaddexp(self.values[key], value)

        def expected(self):
            return sum(key * exp(value) for key, value in self.values.items())

    stats = defaultdict(lambda: defaultdict(ExpectDict))

    for cursor in traversal.depth(circ, unique="id"):
        node = cursor.node
        key = id(node)
        operator = node.data.operator
        prob = boltz_in[key] + boltz_out[key] - boltz_in[id(circ)]

        if operator in ("gain", "loss"):
            info = node.data.args[0]
            stats[info.name][operator].add(info.size, prob)
        elif operator == "split":
            info = node.data.args[0]
            stats[info.name]["out"].add(info.out, prob)
            stats[info.name]["in_left"].add(info.in_left, prob)
            stats[info.name]["in_right"].add(info.in_right, prob)

    for name in stats:
        gain_exp = stats[name]["gain"].expected()
        loss_exp = stats[name]["loss"].expected()
        out_exp = stats[name]["out"].expected()
        in_left_exp = stats[name]["in_left"].expected()
        in_right_exp = stats[name]["in_right"].expected()
        size = len(lca_min[name].node.data.min) + out_exp + in_left_exp + in_right_exp
        print(name)
        print(f"loss={loss_exp:.6f}")
        print(f"gain={gain_exp:.6f}")
        print(f"size={size:.6f}")
        print()


def parse_args():
    """Process the command-line arguments."""
    parser = ArgumentParser(
        description=(
            "Solve small parsimony problems on sets using the InOutParsimony algorithm."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help=(
            "path to the input file describing the problem instance using the NHX "
            "format; use - to read from standard input (default: -)"
        ),
        default="-",
    )

    def sentence_to_description(sentence):
        return sentence[0].lower() + sentence[1:].removesuffix(".")

    subparsers = parser.add_subparsers(title="modes", dest="subcommand", required=True)

    pareto_help = "Compute the set of Pareto-optimal vectors for the given instance."
    subparsers.add_parser(
        "pareto", help=sentence_to_description(pareto_help), description=pareto_help
    )

    sols_help = (
        "Compute information about the set of solutions for the given instance "
        "under a given cost model."
    )
    sols_parser = subparsers.add_parser(
        "solutions",
        help=sentence_to_description(sols_help),
        description=sols_help,
    )
    sols_parser.add_argument(
        "--gain-open",
        "-G",
        type=float,
        metavar="COST",
        default=1,
        help="non-negative cost for each gain event in the solution (default: 1)",
    )
    sols_parser.add_argument(
        "--gain-ext",
        "-E",
        type=float,
        metavar="COST",
        default=0,
        help="non-negative cost for each gained gene in the solution (default: 0)",
    )
    sols_parser.add_argument(
        "--loss-open",
        "-L",
        type=float,
        metavar="COST",
        default=1,
        help="non-negative cost for each loss event in the solution (default: 1)",
    )
    sols_parser.add_argument(
        "--loss-ext",
        "-F",
        type=float,
        metavar="COST",
        default=0,
        help="non-negative cost for each lost gene in the solution (default: 0)",
    )
    sols_parser.add_argument(
        "--temperature",
        "-T",
        type=float,
        metavar="TEMP",
        default=0,
        help=(
            "for `avg` and `show` submodes, select between a distribution including "
            "only optimal solutions (temperature of 0), a distribution including all "
            "solutions (temperature of +inf), or intermediate distributions "
            "(default: 0)"
        ),
    )

    sols_subparsers = sols_parser.add_subparsers(
        title="submodes", dest="solutions", required=True
    )

    avg_help = (
        "Compute the average number of genes in the set at each node of a solution."
    )
    sols_subparsers.add_parser(
        "avg", help=sentence_to_description(avg_help), description=avg_help
    )

    show_help = "Show some or all optimal solutions under the given cost model."
    show_parser = sols_subparsers.add_parser(
        "show", help=sentence_to_description(show_help), description=show_help
    )
    show_parser.add_argument(
        "-c",
        "--count",
        type=float,
        metavar="SIZE",
        default=1,
        help=(
            "number of randomly-sampled solutions to present, "
            "or +inf to present all solutions (default: 1)"
        ),
    )
    show_parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=42,
        help=(
            "seed to use for random generation, in cases where non-deterministic "
            "solution sampling is performed (i.e., when both the temperature and "
            "the sample size are neither 0 nor +inf) (default: 42)"
        ),
    )

    tables_help = (
        "Show the internal dynamic programming tables resulting from running the "
        "algorithm under the minimum cost algebra."
    )
    sols_subparsers.add_parser(
        "tables", help=sentence_to_description(tables_help), description=tables_help
    )

    return parser.parse_args()


def main():
    """Run the InOutParsimony command-line interface."""
    options = parse_args()
    input_file_path = 0 if options.input == "-" else options.input

    with open(input_file_path) as input_file:
        input_tree = read_tree(input_file.read())
        lca_min = compute_lca_min(input_tree)

    if options.subcommand == "pareto":
        do_pareto(lca_min)

    if options.subcommand == "solutions":
        costs = {
            "loss": lambda event: (
                (options.loss_open + options.loss_ext * event.size)
                if event.size > 0
                else 0
            ),
            "gain": lambda event: (
                (options.gain_open + options.gain_ext * event.size)
                if event.size > 0
                else 0
            ),
        }

        if options.solutions == "show":
            do_solutions_show(
                lca_min,
                costs,
                options.temperature,
                Random(options.seed),
                options.count,
            )

        if options.solutions == "tables":
            do_solutions_tables(lca_min, costs)

        if options.solutions == "avg":
            do_solutions_avg(lca_min, costs, options.temperature)
