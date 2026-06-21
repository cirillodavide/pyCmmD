"""Harmonize all layer node identifiers to a common HGNC-based gene index.

Builds a resolver from the HGNC complete set that maps any identifier
(Entrez / UniProt / HGNC id / approved symbol / previous & alias symbols) to the
canonical (approved symbol, hgnc_id). Each layer's edges are then remapped:

  * canonical node label  = HGNC approved symbol
  * source_id / target_id = hgnc_id   (the stable common index)

Resolution priority (most reliable first): the layer's native id -> approved
symbol -> unambiguous previous/alias symbol. Unresolvable nodes are KEPT under
their original symbol with an empty hgnc_id (flagged) and listed in a report.
Edges that collapse to a self-loop after synonym merging are dropped.

Run:  python -m pycmmd.harmonize [--force]
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from . import config
from .utils import download, write_provenance

HGNC_URL = ("https://storage.googleapis.com/public-download-files/hgnc/"
            "tsv/tsv/hgnc_complete_set.txt")

# 0-indexed HGNC complete-set columns.
H_HGNC, H_SYMBOL, H_LOCUS, H_STATUS = 0, 1, 4, 5
H_ALIAS, H_PREV, H_ENTREZ, H_UNIPROT = 8, 10, 18, 25

# native identifier type carried in each layer's source_id/target_id column.
LAYER_NATIVE = {
    "BioGRID": "entrez",
    "ChEMBL": "uniprot",
    "Monarch": "hgnc",
    "Reactome": "entrez",
    "Recon3D": "entrez",
}


class Resolver:
    """Maps identifiers to canonical (symbol, hgnc_id)."""

    def __init__(self, hgnc_path: Path):
        self.by_entrez: dict[str, tuple] = {}
        self.by_uniprot: dict[str, tuple] = {}
        self.by_hgnc: dict[str, tuple] = {}
        self.by_symbol: dict[str, tuple] = {}
        syn_owners: dict[str, set] = defaultdict(set)
        syn_map: dict[str, tuple] = {}
        self.n_genes = 0

        with open(hgnc_path, encoding="utf-8") as fh:
            next(fh)  # header
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) <= H_UNIPROT or p[H_STATUS] != "Approved":
                    continue
                self.n_genes += 1
                symbol, hgnc = p[H_SYMBOL], p[H_HGNC]
                canon = (symbol, hgnc)
                self.by_hgnc[hgnc] = canon
                self.by_symbol[symbol] = canon
                if p[H_ENTREZ]:
                    self.by_entrez[p[H_ENTREZ]] = canon
                for u in p[H_UNIPROT].split("|"):
                    if u:
                        self.by_uniprot[u] = canon
                for s in p[H_ALIAS].split("|") + p[H_PREV].split("|"):
                    if s:
                        syn_owners[s].add(hgnc)
                        syn_map[s] = canon
        # keep only unambiguous synonyms that are not themselves approved symbols
        self.by_synonym = {s: syn_map[s] for s, owners in syn_owners.items()
                           if len(owners) == 1 and s not in self.by_symbol}
        self._native = {"entrez": self.by_entrez,
                        "uniprot": self.by_uniprot, "hgnc": self.by_hgnc}

    def resolve(self, native_type: str, native_id: str, symbol: str):
        """Return (canonical_symbol, hgnc_id, method)."""
        nm = self._native[native_type]
        if native_id and native_id in nm:
            return (*nm[native_id], native_type)
        if symbol in self.by_symbol:
            return (*self.by_symbol[symbol], "symbol")
        if symbol in self.by_synonym:
            return (*self.by_synonym[symbol], "synonym")
        return (symbol, "", "unmapped")


def _harmonize_layer(layer: str, resolver: Resolver) -> dict:
    native = LAYER_NATIVE[layer]
    src = config.PROCESSED_DIR / f"{layer}.edges.tsv"
    out = config.HARMONIZED_DIR / f"{layer}.edges.tsv"
    unmapped_path = config.HARMONIZED_DIR / f"{layer}.unmapped.tsv"

    seen: set = set()
    nodes_in: set = set()
    nodes_out: set = set()
    methods: dict[str, int] = defaultdict(int)
    unmapped: dict[str, str] = {}   # native_id-or-symbol -> original symbol
    edges_in = self_loops = 0

    with open(src, encoding="utf-8") as fh, \
         open(out, "w", newline="", encoding="utf-8") as fo:
        r = csv.reader(fh, delimiter="\t")
        w = csv.writer(fo, delimiter="\t")
        header = next(r)
        w.writerow(config.EDGE_COLUMNS)
        for row in r:
            edges_in += 1
            s_sym, t_sym, lyr, inter, annot, s_id, t_id = row
            nodes_in.update((s_sym, t_sym))
            cs, ch, ms = resolver.resolve(native, s_id, s_sym)
            ct, th, mt = resolver.resolve(native, t_id, t_sym)
            methods[ms] += 1
            methods[mt] += 1
            if ms == "unmapped":
                unmapped[s_id or s_sym] = s_sym
            if mt == "unmapped":
                unmapped[t_id or t_sym] = t_sym
            # order undirected by canonical symbol; keep ids aligned
            a, b, ai, bi = cs, ct, ch, th
            if b < a:
                a, b, ai, bi = ct, cs, th, ch
            if a == b:
                self_loops += 1
                continue
            key = (a, b, inter, annot)
            if key in seen:
                continue
            seen.add(key)
            w.writerow([a, b, lyr, inter, annot, ai, bi])
            nodes_out.update((a, b))

    with open(unmapped_path, "w", newline="", encoding="utf-8") as fu:
        w = csv.writer(fu, delimiter="\t")
        w.writerow(["identifier", "original_symbol"])
        for k, v in sorted(unmapped.items()):
            w.writerow([k, v])

    return {
        "layer": layer,
        "edges_in": edges_in,
        "edges_out": len(seen),
        "self_loops_collapsed": self_loops,
        "nodes_in": len(nodes_in),
        "nodes_out": len(nodes_out),
        "resolve_methods": dict(methods),
        "unmapped_identifiers": len(unmapped),
    }


def build(force: bool = False) -> Path:
    config.ensure_dirs()
    hgnc_path = download(HGNC_URL, config.RAW_DIR / "hgnc_complete_set.txt", force=force)
    print("[parse] building HGNC resolver ...")
    resolver = Resolver(hgnc_path)
    print(f"[parse] {resolver.n_genes:,} approved HGNC genes  "
          f"(entrez={len(resolver.by_entrez):,} uniprot={len(resolver.by_uniprot):,} "
          f"synonyms={len(resolver.by_synonym):,})")

    reports = []
    all_nodes: dict[str, dict] = {}   # symbol -> {hgnc, layers}
    for layer in LAYER_NATIVE:
        if not (config.PROCESSED_DIR / f"{layer}.edges.tsv").exists():
            print(f"[skip] {layer}: no processed edges")
            continue
        rep = _harmonize_layer(layer, resolver)
        reports.append(rep)
        print(f"[{layer}] edges {rep['edges_in']:,} -> {rep['edges_out']:,}  "
              f"nodes {rep['nodes_in']:,} -> {rep['nodes_out']:,}  "
              f"self-loops {rep['self_loops_collapsed']:,}  "
              f"unmapped {rep['unmapped_identifiers']:,}")
        # accumulate node index from harmonized output
        with open(config.HARMONIZED_DIR / f"{layer}.edges.tsv", encoding="utf-8") as fh:
            rr = csv.reader(fh, delimiter="\t")
            next(rr)
            for a, b, lyr, *_rest, ai, bi in (row for row in rr):
                for s, h in ((a, ai), (b, bi)):
                    e = all_nodes.setdefault(s, {"hgnc": h, "layers": set()})
                    e["layers"].add(lyr)
                    if h and not e["hgnc"]:
                        e["hgnc"] = h

    # combined node index
    nodes_path = config.HARMONIZED_DIR / "nodes.tsv"
    with open(nodes_path, "w", newline="", encoding="utf-8") as fn:
        w = csv.writer(fn, delimiter="\t")
        w.writerow(["symbol", "hgnc_id", "harmonized", "n_layers", "layers"])
        for s in sorted(all_nodes):
            e = all_nodes[s]
            w.writerow([s, e["hgnc"], bool(e["hgnc"]), len(e["layers"]),
                        ";".join(sorted(e["layers"]))])

    write_provenance(config.HARMONIZED_DIR / "harmonization_report.json", {
        "hgnc_source": HGNC_URL,
        "hgnc_genes": resolver.n_genes,
        "policy": {"unmapped": "keep_flagged", "locus_scope": "all_approved"},
        "layers": reports,
        "total_unique_nodes": len(all_nodes),
        "harmonized_nodes": sum(1 for e in all_nodes.values() if e["hgnc"]),
    })
    print(f"[ok] wrote harmonized layers + node index to {config.HARMONIZED_DIR}")
    print(f"[ok] common index: {len(all_nodes):,} unique nodes "
          f"({sum(1 for e in all_nodes.values() if e['hgnc']):,} harmonized)")
    return config.HARMONIZED_DIR


def main() -> None:
    ap = argparse.ArgumentParser(description="Harmonize layers to a common HGNC gene index.")
    ap.add_argument("--force", action="store_true", help="re-download HGNC even if cached")
    args = ap.parse_args()
    build(force=args.force)


if __name__ == "__main__":
    main()
