from __future__ import annotations

from functools import lru_cache
from math import log2
from typing import FrozenSet, Iterable, Optional

import numpy as np
from . import utils

Taxon = str
Part = FrozenSet[Taxon]


class SplitInfo:
    ROOT_DUMMY_TAXON = "__ROOT__"
    ROOT_DUMMY_INT = -1

    def information_taxa(self, taxa: Iterable[Taxon], rooted: bool = False):
        """
        Taxon set on which phylogenetic-information calculations should be made.

        For rooted clades, rooted trees are handled by attaching a dummy leaf
        to the root and treating the result as an unrooted split system.
        """
        taxa = utils._to_fset(taxa)
        if not rooted:
            return taxa

        if all(isinstance(t, int) for t in taxa):
            if self.ROOT_DUMMY_INT in taxa:
                raise ValueError(
                    f"Reserved dummy-root integer {self.ROOT_DUMMY_INT!r} already present"
                )
            return frozenset(set(taxa) | {self.ROOT_DUMMY_INT})

        if self.ROOT_DUMMY_TAXON in taxa:
            raise ValueError(
                f"Reserved dummy-root taxon name {self.ROOT_DUMMY_TAXON!r} already present"
            )
        return frozenset(set(taxa) | {self.ROOT_DUMMY_TAXON})


    @staticmethod
    def _mask_complement(mask: int, full_mask: int) -> int:
        return full_mask ^ mask

    @staticmethod
    def _mask_size(mask: int) -> int:
        return mask.bit_count()

    @staticmethod
    def _split_entropy_from_sizes(split_size: int, taxa_size: int) -> float:
        p = split_size / taxa_size
        q = 1.0 - p
        if p == 0.0 or q == 0.0:
            return 0.0
        return -(p * log2(p) + q * log2(q))

    @staticmethod
    def _split_phylogenetic_information_from_sizes(split_size: int, taxa_size: int) -> float:
        """
        h(S) = -log2 PPhy(S), using only split sizes.
        """
        a = split_size
        b = taxa_size - a

        if a < 2 or b < 2:
            return float("inf")

        log_p = (
            utils.log_double_factorial_odd(2 * a - 3)
            + utils.log_double_factorial_odd(2 * b - 3)
            - utils.log_double_factorial_odd(2 * taxa_size - 5)
        )
        return -log_p / np.log(2.0)

    @staticmethod
    def _compatible_unrooted_masks(a_mask: int, b_mask: int, full_mask: int) -> bool:
        """
        Unrooted split compatibility using bitmasks.
        """
        ac = full_mask ^ a_mask
        bc = full_mask ^ b_mask
        return (
            (a_mask & b_mask) == 0
            or (a_mask & bc) == 0
            or (ac & b_mask) == 0
            or (ac & bc) == 0
        )

    @staticmethod
    def _choose_orientation_subset(a_mask: int, b_mask: int, full_mask: int):
        """
        Choose an orientation such that A2 ⊆ A1, if possible.
        Returns (A1, B1, A2, B2) as masks or None.
        """
        ac = full_mask ^ a_mask
        bc = full_mask ^ b_mask

        orientations = [
            (a_mask, ac, b_mask, bc),
            (a_mask, ac, bc, b_mask),
            (ac, a_mask, b_mask, bc),
            (ac, a_mask, bc, b_mask),
        ]

        for A1, B1, A2, B2 in orientations:
            if (A2 & ~A1) == 0:
                return A1, B1, A2, B2

        return None

    def split_entropy(self, split: Part, taxa: Iterable[Taxon]) -> float:
        taxa = utils._to_fset(taxa)
        return self._split_entropy_from_sizes(len(split), len(taxa))

    def split_phylogenetic_information(self, split: Part, taxa: Iterable[Taxon]) -> float:
        taxa = utils._to_fset(taxa)
        return self._split_phylogenetic_information_from_sizes(len(split), len(taxa))

    def clade_entropy(self, clade: Part, taxa: Iterable[Taxon]) -> float:
        taxa = utils._to_fset(taxa)
        return self._split_entropy_from_sizes(len(clade), len(taxa))

    def clade_phylogenetic_information(self, clade: Part, taxa: Iterable[Taxon]) -> float:
        taxa = utils._to_fset(taxa)
        return self._split_phylogenetic_information_from_sizes(len(clade), len(taxa))


    def _compatible_unrooted(self, a: Part, b: Part, taxa: Part) -> bool:
        taxa = utils._to_fset(taxa)
        masks, full_mask, _ = utils._parts_to_masks([a, b], taxa)
        return self._compatible_unrooted_masks(masks[0], masks[1], full_mask)

    @lru_cache(maxsize=utils.CACHE_SIZE_PAIRS)
    def _shared_phylogenetic_information_cached(self, a: Part, b: Part, taxa: Part) -> float:
        """
        Backward-compatible set-based cache wrapper.
        """
        taxa = utils._to_fset(taxa)
        masks, full_mask, _ = utils._parts_to_masks([a, b], taxa)
        return self._shared_phylogenetic_information_masks(masks[0], masks[1], full_mask)

    def _shared_phylogenetic_information_masks(self, a_mask: int, b_mask: int, full_mask: int) -> float:
        """
        Bitmask-native SPI.
        """
        if not self._compatible_unrooted_masks(a_mask, b_mask, full_mask):
            return 0.0

        if a_mask == b_mask:
            return self._split_phylogenetic_information_from_sizes(
                a_mask.bit_count(),
                full_mask.bit_count(),
            )

        chosen = self._choose_orientation_subset(a_mask, b_mask, full_mask)
        if chosen is None:
            return 0.0

        A1, B1, A2, B2 = chosen
        del B2

        n = full_mask.bit_count()
        len_B1 = B1.bit_count()
        len_A2 = A2.bit_count()
        len_A1 = A1.bit_count()

        part1 = 2 * (len_B1 + 1) - 5
        part2 = 2 * (len_A2 + 1) - 5
        part3 = 2 * (len_A1 - len_A2 + 2) - 5
        denom = 2 * n - 5

        log_p12 = (
            utils.log_double_factorial_odd(part1)
            + utils.log_double_factorial_odd(part2)
            + utils.log_double_factorial_odd(part3)
            - utils.log_double_factorial_odd(denom)
        )

        h1 = self._split_phylogenetic_information_from_sizes(a_mask.bit_count(), n)
        h2 = self._split_phylogenetic_information_from_sizes(b_mask.bit_count(), n)
        h12 = -log_p12 / np.log(2.0)

        return max(0.0, h1 + h2 - h12)

    def shared_phylogenetic_information(self, a: Part, b: Part, taxa: Iterable[Taxon]) -> float:
        taxa = utils._to_fset(taxa)
        return self._shared_phylogenetic_information_cached(a, b, taxa)


    def _split_mutual_clustering_information_masks(self, a_mask: int, b_mask: int, full_mask: int) -> float:

        n11 = (a_mask & b_mask).bit_count()
        n10 = (a_mask & (full_mask ^ b_mask)).bit_count()
        n01 = ((full_mask ^ a_mask) & b_mask).bit_count()
        n00 = ((full_mask ^ a_mask) & (full_mask ^ b_mask)).bit_count()

        total = n11 + n10 + n01 + n00
        if total == 0:
            return 0.0

        cells = np.array([[n11, n10], [n01, n00]], dtype=float)
        row = cells.sum(axis=1, keepdims=True)
        col = cells.sum(axis=0, keepdims=True)

        mi = 0.0
        for i in range(2):
            for j in range(2):
                if cells[i, j] == 0:
                    continue
                pij = cells[i, j] / total
                pi = row[i, 0] / total
                pj = col[0, j] / total
                mi += pij * log2(pij / (pi * pj))

        return float(mi)

    def split_mutual_clustering_information(self, a: Part, b: Part, taxa: Iterable[Taxon]) -> float:
        taxa = utils._to_fset(taxa)
        masks, full_mask, _ = utils._parts_to_masks([a, b], taxa)
        return self._split_mutual_clustering_information_masks(masks[0], masks[1], full_mask)


    def _most_informative_compatible_split_mask_unrooted(
        self,
        a_mask: int,
        b_mask: int,
        full_mask: int,
    ) -> Optional[tuple[int, int]]:
        """
        Returns (best_split_mask, candidate_taxa_mask) for MSI.
        """
        ac = full_mask ^ a_mask
        bc = full_mask ^ b_mask

        best = None

        # Candidate 1: (A1∩A2 | B1∩B2)
        x1 = a_mask & b_mask
        y1 = ac & bc
        if x1.bit_count() >= 2 and y1.bit_count() >= 2:
            cand_taxa = x1 | y1
            score = self._split_phylogenetic_information_from_sizes(
                x1.bit_count(), cand_taxa.bit_count()
            )
            best = (score, x1, cand_taxa)

        # Candidate 2: (A1∩B2 | B1∩A2)
        x2 = a_mask & bc
        y2 = ac & b_mask
        if x2.bit_count() >= 2 and y2.bit_count() >= 2:
            cand_taxa = x2 | y2
            score = self._split_phylogenetic_information_from_sizes(
                x2.bit_count(), cand_taxa.bit_count()
            )
            if best is None or score > best[0]:
                best = (score, x2, cand_taxa)

        if best is None:
            return None

        _, best_split, best_taxa = best
        return best_split, best_taxa

    def most_informative_compatible_split_unrooted(
        self,
        a: Part,
        b: Part,
        taxa: Part,
    ) -> Optional[Part]:
        """
        Backward-compatible set-based wrapper.
        """
        taxa = utils._to_fset(taxa)
        masks, full_mask, _ = utils._parts_to_masks([a, b], taxa)
        best = self._most_informative_compatible_split_mask_unrooted(
            masks[0], masks[1], full_mask
        )
        if best is None:
            return None

        best_mask, _ = best
        return frozenset(i for i in range(full_mask.bit_count()) if (best_mask >> i) & 1)

    def most_informative_compatible_clade_rooted(self, a: Part, b: Part) -> Optional[Part]:
        if a <= b:
            return a if len(a) >= 2 else None
        if b <= a:
            return b if len(b) >= 2 else None
        return None