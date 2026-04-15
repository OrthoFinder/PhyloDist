from __future__ import annotations

import os
import re
from typing import FrozenSet, Iterable, List, Optional, Tuple
from ete4 import Tree
from . import utils

Taxon = str
Part = FrozenSet[Taxon]


class TreeStructure:
    """
    Thin wrapper around an ete4.Tree that exposes rooted clades or unrooted splits.

    Accepted input sources:
    - Newick string
    - path to a file containing Newick
    - path to a NEXUS file containing a tree
    - file-like object with .read()
    - ete4.Tree
    """

    def __init__(
        self,
        source: str | Tree,
        is_rooted: Optional[bool] = None,
        parser: int | None = None,
    ):
        self.source = source

        if isinstance(source, Tree):
            self.t = source.copy()
        else:
            text = self._coerce_tree_text(source)
            text = self._clean_tree_text(text)

            if parser is not None:
                self.t = Tree(text, parser=parser)
            else:
                try:
                    self.t = Tree(text, parser=1)
                except Exception:
                    self.t = Tree(text, parser=0)

        self.taxa: FrozenSet[Taxon] = frozenset(self.t.leaf_names())
        self.is_rooted: Optional[bool] = is_rooted

    @classmethod
    def from_newick(
        cls,
        newick: str,
        is_rooted: Optional[bool] = None,
        parser: int | None = None,
    ) -> "TreeStructure":
        return cls(newick, is_rooted=is_rooted, parser=parser)

    @classmethod
    def from_tree(
        cls,
        tree: Tree,
        is_rooted: Optional[bool] = None,
        parser: int | None = None,
    ) -> "TreeStructure":
        return cls(tree, is_rooted=is_rooted, parser=parser)

    @staticmethod
    def _coerce_tree_text(source) -> str:
        """
        Convert input source to text.

        Accepts:
        - file-like object
        - path string
        - raw Newick / NEXUS string
        """
        if hasattr(source, "read"):
            return source.read()

        text = str(source)

        if os.path.exists(text):
            with open(text, "r", encoding="utf-8") as f:
                return f.read()

        return text

    @staticmethod
    def _clean_tree_text(text: str) -> str:
        """
        Normalize tree text into a clean Newick string.

        Supports:
        - plain Newick
        - NEXUS with a tree assignment line
        - optional BEAST-style [&R] prefix
        """
        text = text.strip()

        # NEXUS detection
        if text.lower().startswith("#nexus") or "begin trees" in text.lower():
            match = re.search(r'=\s*(.+?);', text, re.S)
            if not match:
                raise ValueError("Could not find tree in NEXUS content.")
            text = match.group(1).strip()

        # Remove leading rootedness annotation such as [&R]
        text = re.sub(r'^\[&R\]\s*', '', text)

        if not text.endswith(";"):
            text += ";"

        return text

    def __repr__(self) -> str:
        kind = (
            "rooted" if self.is_rooted is True
            else "unrooted" if self.is_rooted is False
            else "unknown-rootedness"
        )
        return f"TreeStructure(n_taxa={len(self.taxa)}, {kind})"

    @property
    def n_taxa(self) -> int:
        return len(self.taxa)

    def get_tree_taxa(self) -> FrozenSet[Taxon]:
        return self.taxa

    def copy(self) -> "TreeStructure":
        return TreeStructure(self.t.copy(), is_rooted=self.is_rooted)

    def _iter_internal_nodes(self):
        for node in self.t.traverse():
            if not node.is_leaf:
                yield node

    def get_rooted_clades(self) -> List[Part]:
        clades: List[Part] = []
        for node in self._iter_internal_nodes():
            leaves = frozenset(node.leaf_names())
            if 1 < len(leaves) < self.n_taxa:
                clades.append(leaves)
        return clades

    def get_unrooted_splits(self) -> List[Part]:
        splits: set[Part] = set()
        for node in self._iter_internal_nodes():
            if node.is_root:
                continue
            leaves = frozenset(node.leaf_names())
            canon = utils._canon_unrooted_split(leaves, self.taxa)
            if canon is not None:
                splits.add(canon)
        return list(splits)

    def normalize_parts(
        self,
        parts: Iterable[Iterable[Taxon]],
        *,
        rooted: bool,
    ) -> List[Part]:
        out: List[Part] = []
        seen: set[Part] = set()

        for p in parts:
            fp = utils._to_fset(p)
            canon = (
                utils._canon_rooted_clade(fp, self.taxa)
                if rooted else
                utils._canon_unrooted_split(fp, self.taxa)
            )
            if canon is not None and canon not in seen:
                seen.add(canon)
                out.append(canon)
        return out

    @staticmethod
    def project_parts_to_shared_taxa(
        parts: Iterable[Part],
        old_taxa: Iterable[Taxon],
        shared_taxa: Iterable[Taxon],
        rooted: bool,
    ) -> List[Part]:
        old_taxa = utils._to_fset(old_taxa)
        shared_taxa = utils._to_fset(shared_taxa)
        del old_taxa

        projected: List[Part] = []
        seen: set[Part] = set()

        for p in parts:
            p = utils._to_fset(p)
            q = p & shared_taxa
            canon = (
                utils._canon_rooted_clade(q, shared_taxa)
                if rooted else
                utils._canon_unrooted_split(q, shared_taxa)
            )
            if canon is not None and canon not in seen:
                seen.add(canon)
                projected.append(canon)

        return projected