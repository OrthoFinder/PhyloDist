from __future__ import annotations

from typing import Optional

from .tree_comparator import TreeComparator
from .tree_structure import TreeStructure
from .utils import ScoreResult


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
    ) -> ScoreResult:
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
            rooted = bool(t1.is_rooted and t2.is_rooted)

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
