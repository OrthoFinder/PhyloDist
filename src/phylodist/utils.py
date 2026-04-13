from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import FrozenSet, Iterable, List, Optional, Tuple, Hashable

import numpy as np
from scipy.special import gammaln

Taxon = Hashable
Part = FrozenSet[Taxon]

CACHE_SIZE_SPLITS = 8192
CACHE_SIZE_PAIRS = 4096


@dataclass
class ScoreResult:
    score: float
    matched_pairs: List[Tuple[Part, Part, float]]
    unmatched1: List[Part]
    unmatched2: List[Part]
    shared_taxa: FrozenSet[Taxon]


DISTANCE_ONLY_METHODS = frozenset({
    "rf",
    "jrf",
    "jrf_bocker",
})

SIMILARITY_ONLY_METHODS = frozenset({
    "nye",
})

DUAL_FORM_METHODS = frozenset({
    "icrf",
    "msi",
    "spi",
    "mci",
})

METHOD_ALIASES = {
    "msid": "msi",
    "pid": "spi",
    "cid": "mci",
}

VALID_RETURN_TYPES = frozenset({"similarity", "distance"})

VALID_METHODS = frozenset({
    "rf", "jrf",
    "jrf_bocker",
    "icrf", "nye",
    "msi", "msid",
    "spi", "pid",
    "mci", "cid",
})

METHOD_SPECS = {
    "rf": {"base": "rf", "fixed_return_type": "distance"},
    "jrf": {"base": "jrf", "fixed_return_type": "distance"},
    "jrf_bocker": {"base": "jrf_bocker", "fixed_return_type": "distance"},

    "icrf": {"base": "icrf", "fixed_return_type": None},
    "nye": {"base": "nye", "fixed_return_type": "similarity"},

    "msi": {"base": "msi", "fixed_return_type": None},
    "msid": {"base": "msi", "fixed_return_type": "distance"},

    "spi": {"base": "spi", "fixed_return_type": None},
    "pid": {"base": "spi", "fixed_return_type": "distance"},

    "mci": {"base": "mci", "fixed_return_type": None},
    "cid": {"base": "mci", "fixed_return_type": "distance"},
}


def _normalize_return_type(return_type: str) -> str:
    return_type = return_type.lower().strip()
    if return_type not in VALID_RETURN_TYPES:
        raise ValueError("return_type must be one of {'similarity', 'distance'}")
    return return_type


def resolve_method_and_return_type(method: str, return_type: str | None = None) -> tuple[str, str]:
    method = method.lower().strip()
    if method not in VALID_METHODS:
        raise ValueError(
            f"Unknown method '{method}'. Expected one of {sorted(VALID_METHODS)}"
        )

    if return_type is not None:
        return_type = _normalize_return_type(return_type)

    spec = METHOD_SPECS[method]
    base = spec["base"]
    fixed = spec["fixed_return_type"]

    if fixed is not None:
        if return_type is None:
            return_type = fixed
        elif return_type != fixed:
            raise ValueError(
                f"Method '{method}' only supports return_type='{fixed}', "
                f"but got return_type='{return_type}'."
            )
        return base, return_type

    if return_type is None:
        return_type = "similarity"

    return base, return_type


def _to_fset(x: Iterable[Taxon]) -> Part:
    return frozenset(x)


def safe_score(x: float) -> float:
    x = float(x)
    return x if np.isfinite(x) else 0.0


@lru_cache(maxsize=CACHE_SIZE_SPLITS)
def log_double_factorial_odd(n: int) -> float:
    """
    Return log(n!!) for odd n >= -1, using gammaln.

    For odd n = 2k - 1:
        (2k - 1)!! = (2k)! / (2^k k!)
    """
    if n <= 1:
        return 0.0
    if n % 2 == 0:
        raise ValueError(f"log_double_factorial_odd expected odd n, got {n}")

    k = (n + 1) // 2
    return gammaln(2 * k + 1) - k * np.log(2.0) - gammaln(k + 1)


@lru_cache(maxsize=CACHE_SIZE_SPLITS)
def _log_pphy_single_cached(split: Part, taxa: Part) -> float:
    a = len(split)
    b = len(taxa) - a

    if a < 2 or b < 2:
        return float("-inf")

    return (
        log_double_factorial_odd(2 * a - 3)
        + log_double_factorial_odd(2 * b - 3)
        - log_double_factorial_odd(2 * len(taxa) - 5)
    )


def _log_pphy_single(split: Part, taxa: Iterable[Taxon]) -> float:
    return _log_pphy_single_cached(split, _to_fset(taxa))


def _canon_unrooted_split(side: Part, taxa: Part) -> Optional[Part]:
    other = taxa - side
    if len(side) < 2 or len(other) < 2:
        return None
    if len(side) < len(other):
        return side
    if len(other) < len(side):
        return other
    return min(side, other, key=lambda s: tuple(sorted(s)))


def _canon_rooted_clade(clade: Part, taxa: Part) -> Optional[Part]:
    if len(clade) < 2 or len(clade) == len(taxa):
        return None
    return clade


def _build_taxon_index(taxa: Part) -> dict[Taxon, int]:
    return {taxon: i for i, taxon in enumerate(sorted(taxa))}


def _part_to_mask(part: Iterable[Taxon], taxon_index: dict[Taxon, int]) -> int:
    part_fs = _to_fset(part)
    missing = [taxon for taxon in part_fs if taxon not in taxon_index]
    if missing:
        raise KeyError(
            "Split contains taxa not present in taxon_index: "
            f"{missing[:10]} (showing up to 10); "
            f"example types: {[type(x).__name__ for x in missing[:10]]}"
        )

    mask = 0
    for taxon in part_fs:
        mask |= 1 << taxon_index[taxon]
    return mask


def _parts_to_masks(
    parts: Iterable[Iterable[Taxon]],
    taxa: Iterable[Taxon],
    taxon_index: dict[Taxon, int] | None = None,
) -> tuple[list[int], int, dict[Taxon, int]]:
    """
    Strict conversion of parts to bitmasks.

    Important:
    - Does NOT silently intersect parts with taxa.
    - Raises KeyError if any taxon in any part is not present in the taxon_index.
    """
    taxa_fs = _to_fset(taxa)
    parts_fs = [_to_fset(p) for p in parts]

    if taxon_index is None:
        taxon_index = _build_taxon_index(taxa_fs)
    else:
        missing_taxa = [t for t in taxa_fs if t not in taxon_index]
        if missing_taxa:
            raise KeyError(
                "Taxa contain entries not present in provided taxon_index: "
                f"{missing_taxa[:10]}"
            )

    masks = [_part_to_mask(p, taxon_index) for p in parts_fs]
    full_mask = (1 << len(taxa_fs)) - 1
    return masks, full_mask, taxon_index


def _mask_to_part(mask: int, taxa_index_to_name: dict[int, Taxon]) -> frozenset[Taxon]:
    out = []
    bit = 0
    while mask:
        if mask & 1:
            out.append(taxa_index_to_name[bit])
        mask >>= 1
        bit += 1
    return frozenset(out)


def _jaccard_from_masks(a: int, b: int) -> float:
    union = (a | b).bit_count()
    if union == 0:
        return 0.0
    inter = (a & b).bit_count()
    return inter / union