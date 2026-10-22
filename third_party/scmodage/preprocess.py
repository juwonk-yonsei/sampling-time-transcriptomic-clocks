"""Python port of the tAge (v1.5.0) preprocessing used by the scaled-diff clocks.

Order per stratum: count filter -> map to mouse Entrez -> RLE -> log10(x+1) ->
per-sample z-scaling -> align to the clock gene list (missing genes = NaN) ->
subtract the per-gene median of the reference samples.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from .paths import TAGE_META


@lru_cache(maxsize=None)
def clock_gene_list() -> tuple[str, ...]:
    genes = (TAGE_META / "Gene_list_rodent_clocks.txt").read_text().split()
    return tuple(str(int(g)) for g in genes)


@lru_cache(maxsize=None)
def _gene_table(species: str) -> pd.DataFrame:
    gt = pd.read_csv(TAGE_META / "metadata" / f"Gene_table_{species}.csv")
    gt = gt[gt["Entrez"].notna()].copy()
    gt["Entrez"] = gt["Entrez"].astype(np.int64).astype(str)
    return gt


@lru_cache(maxsize=None)
def _ortholog_map(species: str) -> pd.Series:
    key = f"Entrez.{species.capitalize()}"
    orth = pd.read_csv(TAGE_META / "metadata" / "Table_of_orthologs.csv", index_col=0)
    orth = orth[orth["Entrez.Mouse"].notna() & orth[key].notna()]
    orth = orth.drop_duplicates(subset=key, keep="first")
    return pd.Series(
        orth["Entrez.Mouse"].astype(np.int64).astype(str).values,
        index=orth[key].astype(np.int64).astype(str).values,
    )


def _detect_id_type(genes: pd.Index, gt: pd.DataFrame) -> tuple[str, pd.Index]:
    genes = genes.astype(str)
    counts = {c: genes.isin(set(gt[c].dropna().astype(str))).sum() for c in ("Ensembl", "Gene.Symbol", "Entrez")}
    stripped = genes.str.replace(r"\.[0-9]+$", "", regex=True)
    n_strip = stripped.isin(set(gt["Ensembl"].dropna())).sum()
    best = max(counts, key=counts.get)
    if n_strip > counts[best]:
        return "Ensembl", stripped
    if counts[best] == 0:
        raise ValueError("no gene identifiers matched the gene table")
    return best, genes


def map_to_mouse_entrez(counts: pd.DataFrame, species: str) -> pd.DataFrame:
    """Rows = gene IDs of `species`, columns = samples. Returns rows = mouse Entrez."""
    gt = _gene_table(species)
    id_type, keys = _detect_id_type(counts.index, gt)
    lut = gt[gt[id_type].notna()].drop_duplicates(subset=id_type, keep="first")
    lut = pd.Series(lut["Entrez"].values, index=lut[id_type].astype(str).values)
    entrez = pd.Series(keys, index=counts.index).map(lut)
    ok = entrez.notna().values
    mat = counts.loc[ok].copy()
    mat.index = entrez[ok].values
    mat = mat.groupby(level=0, sort=True).sum()
    if species == "mouse":
        return mat
    omap = _ortholog_map(species)
    mouse = pd.Series(mat.index, index=mat.index).map(omap)
    ok = mouse.notna().values
    mat = mat.loc[ok]
    mouse = mouse[ok].values
    keep = ~pd.Index(mouse).duplicated(keep="first")
    mat = mat.loc[keep]
    mat.index = mouse[keep]
    return mat


def filter_genes(counts: pd.DataFrame, count_threshold: float = 10, percent_threshold: float = 20) -> pd.DataFrame:
    passing = (counts >= count_threshold).sum(axis=1) >= counts.shape[1] * percent_threshold / 100
    return counts.loc[passing]


def rle_size_factors(counts: pd.DataFrame) -> pd.Series:
    """edgeR::calcNormFactors(method='RLE') times library size (effective library)."""
    x = counts.values.astype(float)
    lib = x.sum(axis=0)
    with np.errstate(divide="ignore"):
        gm = np.exp(np.log(x).mean(axis=1))
    pos = gm > 0
    if pos.sum() == 0:
        raise ValueError("RLE: no gene has a positive count in every sample")
    f = np.median(x[pos] / gm[pos, None], axis=0) / lib
    f = f / np.exp(np.mean(np.log(f)))
    return pd.Series(lib * f, index=counts.columns)


def rle_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    return counts / rle_size_factors(counts) * 1e7


def scale_samples(log_expr: pd.DataFrame) -> pd.DataFrame:
    """R `scale()` on a genes x samples matrix: z-score each sample across genes (ddof=1)."""
    return (log_expr - log_expr.mean(axis=0)) / log_expr.std(axis=0, ddof=1)


def preprocess_stratum(
    counts: pd.DataFrame,
    species: str,
    reference: pd.Series | None = None,
    count_threshold: float = 10,
    percent_threshold: float = 20,
    gene_list: tuple[str, ...] | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Raw counts (genes x samples) of one stratum -> scaled-diff matrix (samples x clock genes).

    `reference` is a boolean Series over samples; its per-gene median is subtracted.
    If None, all samples are the reference (tAge default).
    """
    gene_list = gene_list or clock_gene_list()
    filt = filter_genes(counts, count_threshold, percent_threshold)
    mapped = map_to_mouse_entrez(filt, species)
    log = np.log10(rle_normalize(mapped) + 1)
    scaled = scale_samples(log)
    aligned = scaled.reindex(list(gene_list))
    ref = aligned.columns if reference is None else reference.index[reference.values.astype(bool)]
    if len(ref) == 0:
        ref = aligned.columns
    diff = aligned.sub(aligned[ref].median(axis=1, skipna=True), axis=0)
    info = {
        "n_input_genes": counts.shape[0],
        "n_filtered_genes": filt.shape[0],
        "n_mapped_genes": mapped.shape[0],
        "n_clock_genes_detected": int(aligned.notna().all(axis=1).sum()),
        "n_reference": len(ref),
    }
    return diff.T, info


def preprocess_by(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    split_by: str | None,
    species: str,
    reference_col: str | None = None,
    reference_label=None,
    **kw,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Preprocess each stratum (e.g. cell type) separately and stack the results."""
    strata = [(None, meta.index)] if split_by is None else list(meta.groupby(split_by).groups.items())
    out, infos = [], []
    for level, idx in strata:
        idx = [i for i in counts.columns if i in set(idx)]
        ref = None
        if reference_col is not None:
            ref = meta.loc[idx, reference_col] == reference_label
        diff, info = preprocess_stratum(counts[idx], species, ref, **kw)
        out.append(diff)
        info["stratum"] = level
        info["n_samples"] = len(idx)
        infos.append(info)
    return pd.concat(out).loc[[c for c in counts.columns if c in set().union(*[set(d.index) for d in out])]], pd.DataFrame(infos)
