from __future__ import annotations

from typing import FrozenSet, Sequence, Hashable

import numpy as np
from scipy.optimize import linear_sum_assignment

from .split_info import SplitInfo
from . import utils

Taxon = Hashable
Part = FrozenSet[Taxon]


class Scores(SplitInfo):
    def _effective_scoring_taxa(self, taxa: Part, rooted: bool) -> Part:
        """
        For rooted comparisons, use dummy-root augmented taxa so rooted clades
        are treated as splits on X ∪ {ROOT_DUMMY}, per Smith.
        """
        return self.information_taxa(taxa, rooted=True) if rooted else utils._to_fset(taxa)

    def _solve_weight_matching(
        self,
        W: np.ndarray,
    ) -> tuple[float, list[tuple[int, int, float]]]:
        """
        Safe Hungarian matching on a similarity matrix.
        Any accidental NaN / +/-inf entries are clamped before solving.
        """
        if W.size == 0:
            return 0.0, []

        W = np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0)

        rows, cols = linear_sum_assignment(-W)

        sim = 0.0
        pairs: list[tuple[int, int, float]] = []
        for r, c in zip(rows, cols):
            w = float(W[r, c])
            if w > 0.0:
                sim += w
                pairs.append((r, c, w))

        return sim, pairs

    def score_icrf(
        self,
        parts1: Sequence[Part],
        parts2: Sequence[Part],
        taxa: Part,
        rooted: bool = False,
    ) -> tuple[float, list[tuple[int, int, float]]]:
        """
        Fast exact-match-only ICRF similarity.
        """
        taxa = utils._to_fset(taxa)

        index2 = {p: j for j, p in enumerate(parts2)}
        sim = 0.0
        pairs: list[tuple[int, int, float]] = []

        for i, p in enumerate(parts1):
            j = index2.get(p)
            if j is None:
                continue

            w = (
                self.clade_phylogenetic_information(p, taxa)
                if rooted
                else self.split_phylogenetic_information(p, taxa)
            )
            if np.isfinite(w) and w > 0.0:
                sim += w
                pairs.append((i, j, float(w)))

        return sim, pairs

    def score_nye(
        self,
        parts1: Sequence[Part],
        parts2: Sequence[Part],
        taxa: Part,
        rooted: bool = False,
    ) -> tuple[float, list[tuple[int, int, float]]]:
        """
        Fast Nye similarity over all pairs.

        Assumes parts1, parts2, and taxa are already in the same label space.
        """
        taxa = utils._to_fset(taxa)
        m = len(parts1)
        n = len(parts2)

        if m == 0 or n == 0:
            return 0.0, []

        W = np.empty((m, n), dtype=float)

        if rooted:
            for i, a in enumerate(parts1):
                for j, b in enumerate(parts2):
                    inter = len(a & b)
                    union = len(a | b)
                    W[i, j] = inter / union if union else 0.0
            return self._solve_weight_matching(W)

        masks1, full_mask, taxon_index = utils._parts_to_masks(parts1, taxa)
        masks2, _, _ = utils._parts_to_masks(parts2, taxa, taxon_index=taxon_index)

        for i, a0 in enumerate(masks1):
            a1 = full_mask ^ a0
            for j, b0 in enumerate(masks2):
                b1 = full_mask ^ b0

                s00 = utils._jaccard_from_masks(a0, b0)
                s11 = utils._jaccard_from_masks(a1, b1)
                s01 = utils._jaccard_from_masks(a0, b1)
                s10 = utils._jaccard_from_masks(a1, b0)

                W[i, j] = max(min(s00, s11), min(s01, s10))

        return self._solve_weight_matching(W)

    def score_msi(
        self,
        parts1: Sequence[Part],
        parts2: Sequence[Part],
        taxa: Part,
        rooted: bool = False,
    ) -> tuple[float, list[tuple[int, int, float]]]:
        score_taxa = self._effective_scoring_taxa(taxa, rooted)
        m = len(parts1)
        n = len(parts2)

        if m == 0 or n == 0:
            return 0.0, []

        masks1, full_mask, taxon_index = utils._parts_to_masks(parts1, score_taxa)
        masks2, _, _ = utils._parts_to_masks(parts2, score_taxa, taxon_index=taxon_index)

        W = np.empty((m, n), dtype=float)

        for i, a_mask in enumerate(masks1):
            for j, b_mask in enumerate(masks2):
                best = self._most_informative_compatible_split_mask_unrooted(
                    a_mask,
                    b_mask,
                    full_mask,
                )
                if best is None:
                    W[i, j] = 0.0
                else:
                    best_split_mask, best_taxa_mask = best
                    w = self._split_phylogenetic_information_from_sizes(
                        best_split_mask.bit_count(),
                        best_taxa_mask.bit_count(),
                    )
                    W[i, j] = utils.safe_score(w)

        return self._solve_weight_matching(W)

    def score_mci(
        self,
        parts1: Sequence[Part],
        parts2: Sequence[Part],
        taxa: Part,
        rooted: bool = False,
    ) -> tuple[float, list[tuple[int, int, float]]]:
        score_taxa = self._effective_scoring_taxa(taxa, rooted)
        m = len(parts1)
        n = len(parts2)

        if m == 0 or n == 0:
            return 0.0, []

        masks1, full_mask, taxon_index = utils._parts_to_masks(parts1, score_taxa)
        masks2, _, _ = utils._parts_to_masks(parts2, score_taxa, taxon_index=taxon_index)

        W = np.empty((m, n), dtype=float)

        for i, a_mask in enumerate(masks1):
            for j, b_mask in enumerate(masks2):
                W[i, j] = utils.safe_score(
                    self._split_mutual_clustering_information_masks(
                        a_mask,
                        b_mask,
                        full_mask,
                    )
                )

        return self._solve_weight_matching(W)

    def score_spi(
        self,
        parts1: Sequence[Part],
        parts2: Sequence[Part],
        taxa: Part,
        rooted: bool = False,
    ) -> tuple[float, list[tuple[int, int, float]]]:
        score_taxa = self._effective_scoring_taxa(taxa, rooted)
        m = len(parts1)
        n = len(parts2)

        if m == 0 or n == 0:
            return 0.0, []

        masks1, full_mask, taxon_index = utils._parts_to_masks(parts1, score_taxa)
        masks2, _, _ = utils._parts_to_masks(parts2, score_taxa, taxon_index=taxon_index)

        W = np.empty((m, n), dtype=float)

        for i, a_mask in enumerate(masks1):
            for j, b_mask in enumerate(masks2):
                w = self._shared_phylogenetic_information_masks(
                    a_mask,
                    b_mask,
                    full_mask,
                )
                W[i, j] = utils.safe_score(w)

        return self._solve_weight_matching(W)

    def score_jrf(
        self,
        parts1: Sequence[Part],
        parts2: Sequence[Part],
        taxa: Part,
        k: float = 1.0,
    ) -> tuple[float, list[tuple[int, int, float]]]:
        """
        TreeDist-compatible default JRF (non-arboreal), using bitmasks.

        Similarity is the optimal matching sum of:
            (Nye pair similarity)^k

        Distance is:
            |parts1| + |parts2| - 2 * similarity
        """
        taxa = utils._to_fset(taxa)
        m = len(parts1)
        n = len(parts2)

        if m == 0 or n == 0:
            return float(m + n), []

        masks1, full_mask, taxon_index = utils._parts_to_masks(parts1, taxa)
        masks2, _, _ = utils._parts_to_masks(parts2, taxa, taxon_index=taxon_index)

        W = np.empty((m, n), dtype=float)

        for i, a0 in enumerate(masks1):
            a1 = full_mask ^ a0
            for j, b0 in enumerate(masks2):
                b1 = full_mask ^ b0

                s00 = utils._jaccard_from_masks(a0, b0)
                s11 = utils._jaccard_from_masks(a1, b1)
                s01 = utils._jaccard_from_masks(a0, b1)
                s10 = utils._jaccard_from_masks(a1, b0)

                nye = max(min(s00, s11), min(s01, s10))
                W[i, j] = nye ** k

        sim, pairs = self._solve_weight_matching(W)
        dist = float(m + n - 2.0 * sim)
        return dist, pairs

    def score_jrf_bocker(
        self,
        parts1: Sequence[Part],
        parts2: Sequence[Part],
        taxa: Part,
        k: float = 1.0,
    ) -> tuple[float, list[tuple[int, int, float]]]:
        """
        Böcker et al. non-arboreal JRF using bitmasks.

            δ_k(Y, Y') = 2 - 2 * J(Y, Y')^k
            δ_k(Y, -) = δ_k(-, Y') = 1

        solved by unconstrained bipartite matching.
        """
        taxa = utils._to_fset(taxa)
        m = len(parts1)
        n = len(parts2)

        if m == 0 and n == 0:
            return 0.0, []
        if m == 0:
            return float(n), []
        if n == 0:
            return float(m), []

        masks1, full_mask, taxon_index = utils._parts_to_masks(parts1, taxa)
        masks2, _, _ = utils._parts_to_masks(parts2, taxa, taxon_index=taxon_index)

        size = m + n
        big = 1e9
        C = np.full((size, size), big, dtype=float)

        for i, a0 in enumerate(masks1):
            a1 = full_mask ^ a0
            for j, b0 in enumerate(masks2):
                b1 = full_mask ^ b0

                s00 = utils._jaccard_from_masks(a0, b0)
                s01 = utils._jaccard_from_masks(a0, b1)
                s10 = utils._jaccard_from_masks(a1, b0)
                s11 = utils._jaccard_from_masks(a1, b1)

                jacc = max(s00, s01, s10, s11)
                C[i, j] = 2.0 - 2.0 * (jacc ** k)

        for i in range(m):
            C[i, n + i] = 1.0

        for j in range(n):
            C[m + j, j] = 1.0

        for i in range(m, m + n):
            for j in range(n, n + m):
                C[i, j] = 0.0

        C = np.nan_to_num(C, nan=big, posinf=big, neginf=big)
        rows, cols = linear_sum_assignment(C)

        dist = float(C[rows, cols].sum())
        pairs: list[tuple[int, int, float]] = []
        for r, c in zip(rows, cols):
            if r < m and c < n and np.isfinite(C[r, c]):
                pairs.append((r, c, float(C[r, c])))

        return dist, pairs