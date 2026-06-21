"""Reusable gene identifier mappings from the official NCBI dictionary.

Downloads NCBI ``Homo_sapiens.gene_info`` (the authoritative Entrez<->Symbol
source, also carrying synonyms) and exposes simple lookup dicts. This is the
backbone the cross-layer harmonization step will build on.
"""

from __future__ import annotations

import gzip

from . import config
from .utils import download

GENE_INFO_URL = ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
                 "Mammalia/Homo_sapiens.gene_info.gz")

# gene_info columns (0-indexed): 1=GeneID(Entrez) 2=Symbol 4=Synonyms
_COL_GENEID = 1
_COL_SYMBOL = 2
_COL_SYNONYMS = 4


def _rows(force: bool = False):
    config.ensure_dirs()
    path = download(GENE_INFO_URL, config.RAW_DIR / "Homo_sapiens.gene_info.gz", force=force)
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        next(fh)  # header
        for line in fh:
            yield line.rstrip("\n").split("\t")


def entrez_to_symbol(force: bool = False) -> dict[str, str]:
    """Return {Entrez GeneID (str) -> official Symbol}."""
    return {p[_COL_GENEID]: p[_COL_SYMBOL] for p in _rows(force)
            if len(p) > _COL_SYMBOL and p[_COL_SYMBOL] != "-"}
