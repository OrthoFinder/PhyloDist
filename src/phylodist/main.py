from __future__ import annotations

from typing import Optional
import argparse
import json

from .tree_comparator import TreeComparator
from .tree_structure import TreeStructure
from . import utils


class PhyloDist:
    """
    Public interface for comparing phylogenetic trees.

    Trees must be provided as TreeStructure objects.
    """

    def __init__(
        self,
        tree1: Optional[TreeStructure] = None,
        tree2: Optional[TreeStructure] = None,
    ):
        self.tree1 = tree1
        self.tree2 = tree2

    def set_trees(
        self,
        tree1: TreeStructure,
        tree2: TreeStructure,
    ) -> "PhyloDist":
        """
        Replace stored trees.
        """
        self.tree1 = tree1
        self.tree2 = tree2
        return self
    
    @classmethod
    def from_newick(
        cls,
        tree1: str,
        tree2: str,
        is_rooted: bool | None = None,
        parser: int | None = None,
    ):
        t1 = TreeStructure(tree1, is_rooted=is_rooted, parser=parser)
        t2 = TreeStructure(tree2, is_rooted=is_rooted, parser=parser)
        return cls(t1, t2)

    def compare(
        self,
        method: str = "msi",
        rooted: Optional[bool] = None,
        return_type: str = "distance",
        normalise: bool = False,
        k: float = 1.0,
        tree1: Optional[TreeStructure] = None,
        tree2: Optional[TreeStructure] = None,
    ) -> utils.ScoreResult:
        """
        Compare two trees.

        Parameters
        ----------
        method
            Comparison method name.

        rooted
            Comparison mode.
            If None, inferred from stored tree metadata.

        return_type
            "similarity" or "distance".

        normalise
            Whether to scale the score to [0,1].

        k
            Exponent parameter used by JRF.

        tree1, tree2
            Optional trees to override stored ones.
        """

        t1 = tree1 if tree1 is not None else self.tree1
        t2 = tree2 if tree2 is not None else self.tree2

        if t1 is None or t2 is None:
            raise ValueError(
                "Two TreeStructure objects are required."
            )

        if not isinstance(t1, TreeStructure):
            raise TypeError("tree1 must be TreeStructure")

        if not isinstance(t2, TreeStructure):
            raise TypeError("tree2 must be TreeStructure")

        if rooted is None:
            if t1.is_rooted is None or t2.is_rooted is None:
                raise ValueError(
                    "Tree rootedness is unknown. "
                    "Specify rooted=True or rooted=False."
                )

            if t1.is_rooted != t2.is_rooted:
                raise ValueError(
                    "Tree rootedness mismatch: "
                    f"tree1 is rooted={t1.is_rooted}, "
                    f"tree2 is rooted={t2.is_rooted}. "
                    "Specify rooted=True or rooted=False."
                )
            
            rooted = t1.is_rooted

        comparator = TreeComparator.from_tree_structures(
            t1,
            t2,
            rooted=rooted,
        )

        return comparator.compare_tree_pair(
            method=method,
            rooted=rooted,
            return_type=return_type,
            normalise=normalise,
            k=k,
        )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Compare phylogenetic trees using PhyloDist."
    )

    parser.add_argument(
        "tree1",
        help="First tree (Newick string or file path)"
    )

    parser.add_argument(
        "tree2",
        help="Second tree (Newick string or file path)"
    )

    parser.add_argument(
        "--method",
        default="rf",
        choices=[
            "rf",
            "jrf",
            "icrf",
            "nye",
            "msi",
            "msid",
            "spi",
            "pid",
            "mci",
            "cid",
        ],
        help="Comparison method"
    )

    parser.add_argument(
        "--normalise",
        action="store_true",
        help="Return normalized score"
    )

    parser.add_argument(
        "--parser",
        type=int,
        default=1,
        help="ETE parser format (default=1)"
    )

    parser.add_argument(
        "--k",
        type=float,
        default=1.0,
        help="Exponent parameter for JRF"
    )

    parser.add_argument(
        "--return-type",
        default="distance",
        choices=["distance", "similarity"],
        help="Return similarity or distance"
    )

    # Tree 1 rootedness metadata
    tree1_root_group = parser.add_mutually_exclusive_group()
    tree1_root_group.add_argument(
        "--tree1-rooted",
        dest="tree1_rooted",
        action="store_true",
        help="Mark tree1 as rooted"
    )
    tree1_root_group.add_argument(
        "--tree1-unrooted",
        dest="tree1_rooted",
        action="store_false",
        help="Mark tree1 as unrooted"
    )
    parser.set_defaults(tree1_rooted=None)

    # Tree 2 rootedness metadata
    tree2_root_group = parser.add_mutually_exclusive_group()
    tree2_root_group.add_argument(
        "--tree2-rooted",
        dest="tree2_rooted",
        action="store_true",
        help="Mark tree2 as rooted"
    )
    tree2_root_group.add_argument(
        "--tree2-unrooted",
        dest="tree2_rooted",
        action="store_false",
        help="Mark tree2 as unrooted"
    )
    parser.set_defaults(tree2_rooted=None)

    # Comparison mode
    compare_root_group = parser.add_mutually_exclusive_group()
    compare_root_group.add_argument(
        "--compare-rooted",
        dest="compare_rooted",
        action="store_true",
        help="Force rooted comparison"
    )
    compare_root_group.add_argument(
        "--compare-unrooted",
        dest="compare_rooted",
        action="store_false",
        help="Force unrooted comparison"
    )
    parser.set_defaults(compare_rooted=None)

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional diagnostic information"
    )

    parser.add_argument(
        "--print-matching",
        action="store_true",
        help="Print matched split/clade pairs"
    )

    parser.add_argument(
        "--json",
        action="store_true",
        help="Output full result as JSON"
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    t1 = TreeStructure(
        args.tree1,
        is_rooted=args.tree1_rooted,
        parser=args.parser,
    )

    t2 = TreeStructure(
        args.tree2,
        is_rooted=args.tree2_rooted,
        parser=args.parser,
    )

    grf = PhyloDist(t1, t2)

    result = grf.compare(
        method=args.method,
        rooted=args.compare_rooted,
        normalise=args.normalise,
        return_type=args.return_type,
        k=args.k,
    )

    if args.json:
        payload = {
            "tree1": args.tree1,
            "tree2": args.tree2,
            "method": args.method,
            "return_type": args.return_type,
            "normalise": args.normalise,
            "k": args.k,
            "tree1_rooted": args.tree1_rooted,
            "tree2_rooted": args.tree2_rooted,
            "compare_rooted": args.compare_rooted,
            "result": utils.serialize_result(result),
        }
        print(json.dumps(payload, indent=2))
        return

    print(result.score)

    if args.verbose:
        print(f"method: {args.method}")
        print(f"return_type: {args.return_type}")
        print(f"normalise: {args.normalise}")
        print(f"k: {args.k}")
        print(f"tree1_rooted: {args.tree1_rooted}")
        print(f"tree2_rooted: {args.tree2_rooted}")
        print(f"compare_rooted: {args.compare_rooted}")
        print(f"shared_taxa: {len(result.shared_taxa)}")
        print(f"matched_pairs: {len(result.matched_pairs)}")
        print(f"unmatched1: {len(result.unmatched1)}")
        print(f"unmatched2: {len(result.unmatched2)}")

    if args.print_matching:
        print("\nMatched pairs:")
        for a, b, w in result.matched_pairs:
            print(
                f"tree1={utils._part_to_jsonable(a)} "
                f"<-> tree2={utils._part_to_jsonable(b)} "
                f"weight={w}"
            )

if __name__ == "__main__":
    main()