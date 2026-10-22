"""Expression -> tAge scaled-diff (samples x clock genes in mouse Entrez space).

ARCHS4 symbols are uppercase. `harmonize_mouse_symbols` matches them to the tAge mouse
gene table without regard to case before the table's own exact-symbol lookup.

RNA-seq counts follow the tAge v1.5.0 order (count filter, RLE, log10, per-sample scaling,
subtraction of the reference median). Microarray data enter on a log scale (one row per gene)
and skip the count-specific steps; per-sample scaling makes the log base irrelevant.
Both return the per-gene mean log expression used to match random genes by expression level.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scmodage import preprocess as pp


def _finish(log: pd.DataFrame, reference: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    scaled = pp.scale_samples(log)
    aligned = scaled.reindex(list(pp.clock_gene_list()))
    ref = reference.index[reference.values.astype(bool)]
    diff = aligned.sub(aligned[ref].median(axis=1, skipna=True), axis=0).T
    level = log.mean(axis=1).reindex(aligned.index)
    return diff, level


def harmonize_mouse_symbols(counts: pd.DataFrame) -> pd.DataFrame:
    """Rename an uppercase symbol index to the tAge gene-table spelling and sum duplicates."""
    gt = pp._gene_table("mouse")
    symbol = gt["Gene.Symbol"].astype(str)
    canon = pd.Series(symbol.values, index=symbol.str.upper()).groupby(level=0).first()
    mapped = pd.Series(pd.Index(counts.index.astype(str)).str.upper().map(canon).to_numpy(), index=counts.index)
    keep = mapped.notna()
    out = counts.loc[keep.to_numpy()].copy()
    out.index = pd.Index(mapped[keep].to_numpy())
    return out.groupby(level=0, sort=False).sum()


def counts_to_diff(counts: pd.DataFrame, species: str, reference: pd.Series,
                   count_threshold: float = 10, percent_threshold: float = 20) -> tuple[pd.DataFrame, pd.Series]:
    """counts: genes (Entrez, Ensembl or symbol of `species`) x samples."""
    filt = pp.filter_genes(counts, count_threshold, percent_threshold)
    mapped = pp.map_to_mouse_entrez(filt, species)
    return _finish(np.log10(pp.rle_normalize(mapped) + 1), reference.loc[counts.columns])


def array_to_diff(log_expr: pd.DataFrame, species: str, reference: pd.Series,
                  min_quantile: float = 0.25) -> tuple[pd.DataFrame, pd.Series]:
    """log_expr: genes (Entrez of `species`, one row per gene) x samples, log scale.
    Genes with mean expression below the `min_quantile` quantile are dropped."""
    mean = log_expr.mean(axis=1)
    m = log_expr.loc[mean >= mean.quantile(min_quantile)]
    if species != "mouse":
        mouse = pd.Series(m.index.astype(str), index=m.index).map(pp._ortholog_map(species))
        m = m.loc[mouse.notna().values]
        m.index = mouse.dropna().values
        m = m[~m.index.duplicated(keep="first")]
    return _finish(m, reference.loc[log_expr.columns])
