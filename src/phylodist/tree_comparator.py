from __future__ import annotations

from typing import FrozenSet, Iterable, Optional, Sequence, Hashable

from .tree_structure import TreeStructure
from .utils import ScoreResult
from . import utils
from .scores_funcs import Scores


Taxon = Hashable
Part = FrozenSet[Taxon]
IntTaxon = int
IntPart = FrozenSet[IntTaxon]


class TreeComparator(Scores):
    def __init__(
        self,
        ts1: Iterable[Iterable[Taxon]],
        ts2: Iterable[Iterable[Taxon]],
        taxa1: Iterable[Taxon],
        taxa2: Iterable[Taxon],
    ):
        self.ts1 = [utils._to_fset(p) for p in ts1]
        self.ts2 = [utils._to_fset(p) for p in ts2]
        self.taxa1 = utils._to_fset(taxa1)
        self.taxa2 = utils._to_fset(taxa2)

    @classmethod
    def from_tree_structures(
        cls,
        tree1: TreeStructure,
        tree2: TreeStructure,
        rooted: bool,
    ) -> "TreeComparator":
        if rooted:
            ts1 = tree1.normalize_parts(tree1.get_rooted_clades(), rooted=True)
            ts2 = tree2.normalize_parts(tree2.get_rooted_clades(), rooted=True)
        else:
            ts1 = tree1.normalize_parts(tree1.get_unrooted_splits(), rooted=False)
            ts2 = tree2.normalize_parts(tree2.get_unrooted_splits(), rooted=False)

        return cls(
            ts1=ts1,
            ts2=ts2,
            taxa1=tree1.get_tree_taxa(),
            taxa2=tree2.get_tree_taxa(),
        )

    def shared_taxa(self) -> FrozenSet[Taxon]:
        return self.taxa1 & self.taxa2

    @staticmethod
    def _build_taxon_codec(
        taxa: Iterable[Taxon],
    ) -> tuple[dict[Taxon, int], dict[int, Taxon], FrozenSet[int]]:
        taxa_fs = utils._to_fset(taxa)
        ordered = sorted(taxa_fs)
        taxon_to_int = {taxon: i for i, taxon in enumerate(ordered)}
        int_to_taxon = {i: taxon for taxon, i in taxon_to_int.items()}
        int_taxa = frozenset(range(len(ordered)))
        return taxon_to_int, int_to_taxon, int_taxa

    @staticmethod
    def _encode_part(
        part: Part,
        taxon_to_int: dict[Taxon, int],
    ) -> IntPart:
        return frozenset(taxon_to_int[t] for t in part)

    @classmethod
    def _encode_parts(
        cls,
        parts: Sequence[Part],
        taxon_to_int: dict[Taxon, int],
    ) -> list[IntPart]:
        return [cls._encode_part(p, taxon_to_int) for p in parts]

    @staticmethod
    def _decode_part(
        part: IntPart,
        int_to_taxon: dict[int, Taxon],
    ) -> Part:
        return frozenset(int_to_taxon[i] for i in part)

    @classmethod
    def _decode_parts(
        cls,
        parts: Sequence[IntPart],
        int_to_taxon: dict[int, Taxon],
    ) -> list[Part]:
        return [cls._decode_part(p, int_to_taxon) for p in parts]

    def similarity_to_distance(
        self,
        similarity: float,
        parts1: Sequence[IntPart],
        parts2: Sequence[IntPart],
        taxa: Iterable[IntTaxon],
        rooted: bool,
        mode: str,
    ) -> float:
        """
        Convert similarity to distance.

        TreeDist-style conversion:
            distance = info(tree1) + info(tree2) - 2 * similarity
        """
        taxa = utils._to_fset(taxa)

        if mode == "mci":
            if rooted:
                total1 = sum(self.clade_entropy(s, taxa) for s in parts1)
                total2 = sum(self.clade_entropy(s, taxa) for s in parts2)
            else:
                total1 = sum(self.split_entropy(s, taxa) for s in parts1)
                total2 = sum(self.split_entropy(s, taxa) for s in parts2)

        elif mode in {"msi", "spi", "icrf"}:
            if rooted:
                total1 = sum(self.clade_phylogenetic_information(s, taxa) for s in parts1)
                total2 = sum(self.clade_phylogenetic_information(s, taxa) for s in parts2)
            else:
                total1 = sum(self.split_phylogenetic_information(s, taxa) for s in parts1)
                total2 = sum(self.split_phylogenetic_information(s, taxa) for s in parts2)
        else:
            raise ValueError(f"Unknown mode {mode}")

        return (total1 + total2) - 2.0 * similarity

    def normalise_score(
        self,
        score: float,
        parts1: Sequence[IntPart],
        parts2: Sequence[IntPart],
        taxa: Iterable[IntTaxon],
        method: str,
        return_type: str,
        rooted: bool = False,
    ) -> float:
        taxa = utils._to_fset(taxa)

        if method in {"rf", "jrf", "jrf_bocker"}:
            denom = len(parts1) + len(parts2)
            return score / denom if denom else 0.0

        if method == "nye":
            denom = min(len(parts1), len(parts2))
            return score / denom if denom else 0.0

        if method == "mci":
            if rooted:
                total1 = sum(self.clade_entropy(s, taxa) for s in parts1)
                total2 = sum(self.clade_entropy(s, taxa) for s in parts2)
            else:
                total1 = sum(self.split_entropy(s, taxa) for s in parts1)
                total2 = sum(self.split_entropy(s, taxa) for s in parts2)

            hmax = 0.5 * (total1 + total2)
            if hmax == 0:
                return 0.0

            if return_type == "similarity":
                return score / hmax
            return score / (total1 + total2)

        if method in {"icrf", "msi", "spi"}:
            if rooted:
                total1 = sum(self.clade_phylogenetic_information(s, taxa) for s in parts1)
                total2 = sum(self.clade_phylogenetic_information(s, taxa) for s in parts2)
            else:
                total1 = sum(self.split_phylogenetic_information(s, taxa) for s in parts1)
                total2 = sum(self.split_phylogenetic_information(s, taxa) for s in parts2)

            hmax = 0.5 * (total1 + total2)
            if hmax == 0:
                return 0.0

            if return_type == "similarity":
                return score / hmax
            return score / (total1 + total2)

        return score

    def compare_tree_pair(
        self,
        p1_normalized: Optional[Iterable[Iterable[Taxon]]] = None,
        p2_normalized: Optional[Iterable[Iterable[Taxon]]] = None,
        method: str = "msi",
        rooted: bool = False,
        k: float = 1.0,
        return_type: str = "distance",
        normalise: bool = False,
    ) -> ScoreResult:
        method, return_type = utils.resolve_method_and_return_type(method, return_type)

        shared = self.shared_taxa()
        if len(shared) < 4:
            return ScoreResult(
                score=float("nan"),
                matched_pairs=[],
                unmatched1=[],
                unmatched2=[],
                shared_taxa=frozenset(),
            )

        p1_input = self.ts1 if p1_normalized is None else [utils._to_fset(p) for p in p1_normalized]
        p2_input = self.ts2 if p2_normalized is None else [utils._to_fset(p) for p in p2_normalized]

        # 1) project both trees to the shared leaf set using original labels
        p1_str = TreeStructure.project_parts_to_shared_taxa(
            p1_input,
            self.taxa1,
            shared,
            rooted=rooted,
        )
        p2_str = TreeStructure.project_parts_to_shared_taxa(
            p2_input,
            self.taxa2,
            shared,
            rooted=rooted,
        )

        taxon_to_int, int_to_taxon, shared_int = self._build_taxon_codec(shared)
        p1 = self._encode_parts(p1_str, taxon_to_int)
        p2 = self._encode_parts(p2_str, taxon_to_int)

        if method == "rf":
            set1 = set(p1)
            set2 = set(p2)
            exact = sorted(set1 & set2, key=lambda s: (len(s), tuple(sorted(s))))
            rf = len(set1 ^ set2)

            score = float(rf) if return_type == "distance" else float(len(exact))

            result = ScoreResult(
                score=score,
                matched_pairs=[
                    (
                        self._decode_part(s, int_to_taxon),
                        self._decode_part(s, int_to_taxon),
                        1.0,
                    )
                    for s in exact
                ],
                unmatched1=[self._decode_part(s, int_to_taxon) for s in p1 if s not in set2],
                unmatched2=[self._decode_part(s, int_to_taxon) for s in p2 if s not in set1],
                shared_taxa=shared,
            )

        elif method == "jrf":
            dist, matched_idx = self.score_jrf(p1, p2, shared_int, k=k)

            matched_i = {i for i, _, _ in matched_idx}
            matched_j = {j for _, j, _ in matched_idx}

            result = ScoreResult(
                score=dist,
                matched_pairs=[(p1_str[i], p2_str[j], w) for i, j, w in matched_idx],
                unmatched1=[p1_str[i] for i in range(len(p1)) if i not in matched_i],
                unmatched2=[p2_str[j] for j in range(len(p2)) if j not in matched_j],
                shared_taxa=shared,
            )

        elif method == "jrf_bocker":
            dist, matched_idx = self.score_jrf_bocker(p1, p2, shared_int, k=k)

            matched_i = {i for i, _, _ in matched_idx}
            matched_j = {j for _, j, _ in matched_idx}

            result = ScoreResult(
                score=dist,
                matched_pairs=[(p1_str[i], p2_str[j], w) for i, j, w in matched_idx],
                unmatched1=[p1_str[i] for i in range(len(p1)) if i not in matched_i],
                unmatched2=[p2_str[j] for j in range(len(p2)) if j not in matched_j],
                shared_taxa=shared,
            )

        elif method == "nye":
            sim, matched_idx = self.score_nye(p1, p2, shared_int, rooted=rooted)

            matched_i = {i for i, _, _ in matched_idx}
            matched_j = {j for _, j, _ in matched_idx}

            result = ScoreResult(
                score=sim,
                matched_pairs=[(p1_str[i], p2_str[j], w) for i, j, w in matched_idx],
                unmatched1=[p1_str[i] for i in range(len(p1)) if i not in matched_i],
                unmatched2=[p2_str[j] for j in range(len(p2)) if j not in matched_j],
                shared_taxa=shared,
            )

        elif method == "icrf":
            sim, matched_idx = self.score_icrf(p1, p2, shared_int, rooted=rooted)

            matched_i = {i for i, _, _ in matched_idx}
            matched_j = {j for _, j, _ in matched_idx}

            score = (
                self.similarity_to_distance(
                    similarity=sim,
                    parts1=p1,
                    parts2=p2,
                    taxa=shared_int,
                    rooted=rooted,
                    mode=method,
                )
                if return_type == "distance"
                else sim
            )

            result = ScoreResult(
                score=score,
                matched_pairs=[(p1_str[i], p2_str[j], w) for i, j, w in matched_idx],
                unmatched1=[p1_str[i] for i in range(len(p1)) if i not in matched_i],
                unmatched2=[p2_str[j] for j in range(len(p2)) if j not in matched_j],
                shared_taxa=shared,
            )

        elif method == "mci":
            sim, matched_idx = self.score_mci(p1, p2, shared_int, rooted=rooted)

            matched_i = {i for i, _, _ in matched_idx}
            matched_j = {j for _, j, _ in matched_idx}

            score = (
                self.similarity_to_distance(
                    similarity=sim,
                    parts1=p1,
                    parts2=p2,
                    taxa=shared_int,
                    rooted=rooted,
                    mode=method,
                )
                if return_type == "distance"
                else sim
            )

            result = ScoreResult(
                score=score,
                matched_pairs=[(p1_str[i], p2_str[j], w) for i, j, w in matched_idx],
                unmatched1=[p1_str[i] for i in range(len(p1)) if i not in matched_i],
                unmatched2=[p2_str[j] for j in range(len(p2)) if j not in matched_j],
                shared_taxa=shared,
            )

        elif method == "spi":
            sim, matched_idx = self.score_spi(p1, p2, shared_int, rooted=rooted)

            matched_i = {i for i, _, _ in matched_idx}
            matched_j = {j for _, j, _ in matched_idx}

            score = (
                self.similarity_to_distance(
                    similarity=sim,
                    parts1=p1,
                    parts2=p2,
                    taxa=shared_int,
                    rooted=rooted,
                    mode=method,
                )
                if return_type == "distance"
                else sim
            )

            result = ScoreResult(
                score=score,
                matched_pairs=[(p1_str[i], p2_str[j], w) for i, j, w in matched_idx],
                unmatched1=[p1_str[i] for i in range(len(p1)) if i not in matched_i],
                unmatched2=[p2_str[j] for j in range(len(p2)) if j not in matched_j],
                shared_taxa=shared,
            )

        elif method == "msi":
            sim, matched_idx = self.score_msi(p1, p2, shared_int, rooted=rooted)

            matched_i = {i for i, _, _ in matched_idx}
            matched_j = {j for _, j, _ in matched_idx}

            score = (
                self.similarity_to_distance(
                    similarity=sim,
                    parts1=p1,
                    parts2=p2,
                    taxa=shared_int,
                    rooted=rooted,
                    mode=method,
                )
                if return_type == "distance"
                else sim
            )

            result = ScoreResult(
                score=score,
                matched_pairs=[(p1_str[i], p2_str[j], w) for i, j, w in matched_idx],
                unmatched1=[p1_str[i] for i in range(len(p1)) if i not in matched_i],
                unmatched2=[p2_str[j] for j in range(len(p2)) if j not in matched_j],
                shared_taxa=shared,
            )

        else:
            raise ValueError(f"Unsupported method: {method}")

        if normalise:
            result.score = self.normalise_score(
                result.score,
                p1,
                p2,
                shared_int,
                method,
                return_type,
                rooted=rooted,
            )

        return result