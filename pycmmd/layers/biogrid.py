"""Layer 1 -- BioGRID molecular interactions.

Downloads the current BioGRID release (tab3 format), keeps human-human
interactions, and writes a standardized edge list.

Edge semantics for this layer:
  * interaction = Experimental System Type  -> "physical" / "genetic"
  * annotation  = Experimental System       -> the specific assay
Nodes are Official Gene Symbols; the stable Entrez Gene ID is retained in the
``source_id`` / ``target_id`` columns to anchor later harmonization.
A companion detail file aggregates, per edge, the supporting BioGRID records:
throughput, record count, publications (PubMed) and BioGRID interaction ids.

Run:  python -m pycmmd.layers.biogrid [--force]
"""

from __future__ import annotations

import argparse
import csv
import io
import zipfile
from pathlib import Path

from .. import config
from ..utils import download, write_provenance

LAYER = "BioGRID"
URL = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.tab3.zip"

# 0-indexed column positions in the BioGRID tab3 format.
COL_BIOGRID_ID = 0         # BioGRID interaction id (detail)
COL_ENTREZ_A = 1
COL_ENTREZ_B = 2
COL_SYMBOL_A = 7
COL_SYMBOL_B = 8
COL_EXP_SYSTEM = 11        # annotation (assay, e.g. "Two-hybrid")
COL_EXP_SYSTEM_TYPE = 12   # interaction ("physical" / "genetic")
COL_PUBLICATION = 14       # publication source, e.g. "PUBMED:12345" (detail)
COL_ORGANISM_A = 15
COL_ORGANISM_B = 16
COL_THROUGHPUT = 17        # Low/High Throughput (detail)

EMPTY = {"-", "", None}


def _open_tab3(zip_path: Path):
    """Yield the inner tab3 .txt member of the BioGRID zip and its version."""
    zf = zipfile.ZipFile(zip_path)
    member = next(n for n in zf.namelist() if n.endswith(".txt"))
    # filename looks like BIOGRID-ALL-4.4.246.tab3.txt -> grab the version token
    version = member.replace("BIOGRID-ALL-", "").split(".tab3")[0]
    stream = io.TextIOWrapper(zf.open(member), encoding="utf-8", newline="")
    return stream, version


def build(force: bool = False, drop_self_loops: bool = True) -> Path:
    config.ensure_dirs()
    zip_path = download(URL, config.RAW_DIR / "BIOGRID-ALL-LATEST.tab3.zip", force=force)

    out_path = config.PROCESSED_DIR / f"{LAYER}.edges.tsv"
    detail_path = config.PROCESSED_DIR / f"{LAYER}.edges.detail.tsv"
    stream, version = _open_tab3(zip_path)
    print(f"[parse] BioGRID version {version}")

    # edge key -> aggregated evidence across the BioGRID records backing it
    edges: dict[tuple, dict] = {}
    rows_total = self_loops = 0

    with stream as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # skip the BioGRID header line
        for row in reader:
            rows_total += 1
            if len(row) <= COL_THROUGHPUT:
                continue
            # human-human only
            if row[COL_ORGANISM_A] != config.HUMAN_TAXON or row[COL_ORGANISM_B] != config.HUMAN_TAXON:
                continue
            sym_a, sym_b = row[COL_SYMBOL_A].strip(), row[COL_SYMBOL_B].strip()
            if sym_a in EMPTY or sym_b in EMPTY:
                continue
            if drop_self_loops and sym_a == sym_b:
                self_loops += 1
                continue

            interaction = row[COL_EXP_SYSTEM_TYPE].strip()  # physical / genetic
            annotation = row[COL_EXP_SYSTEM].strip()        # assay
            ent_a, ent_b = row[COL_ENTREZ_A].strip(), row[COL_ENTREZ_B].strip()

            # undirected dedup, distinct per (interaction, annotation)
            a, b, ia, ib = (sym_a, sym_b, ent_a, ent_b)
            if b < a:
                a, b, ia, ib = sym_b, sym_a, ent_b, ent_a
            key = (a, b, interaction, annotation)

            ev = edges.get(key)
            if ev is None:
                ev = edges[key] = {"ia": ia, "ib": ib, "n": 0,
                                   "pubs": set(), "thr": set(), "bgid": set()}
            ev["n"] += 1
            if row[COL_PUBLICATION] not in EMPTY:
                ev["pubs"].add(row[COL_PUBLICATION].strip())
            if row[COL_THROUGHPUT] not in EMPTY:
                ev["thr"].add(row[COL_THROUGHPUT].strip())
            if row[COL_BIOGRID_ID] not in EMPTY:
                ev["bgid"].add(row[COL_BIOGRID_ID].strip())

    type_counts: dict[str, int] = {}
    nodes: set = set()
    with open(out_path, "w", newline="", encoding="utf-8") as out, \
         open(detail_path, "w", newline="", encoding="utf-8") as det:
        w = csv.writer(out, delimiter="\t")
        dw = csv.writer(det, delimiter="\t")
        w.writerow(config.EDGE_COLUMNS)
        dw.writerow(["source", "target", "experimental_system_type", "experimental_system",
                     "throughput", "n_records", "publications", "biogrid_interaction_ids"])
        for (a, b, interaction, annotation), ev in edges.items():
            w.writerow([a, b, LAYER, interaction, annotation, ev["ia"], ev["ib"]])
            dw.writerow([a, b, interaction, annotation,
                         ";".join(sorted(ev["thr"])), ev["n"],
                         ";".join(sorted(ev["pubs"])), ";".join(sorted(ev["bgid"]))])
            type_counts[interaction] = type_counts.get(interaction, 0) + 1
            nodes.update((a, b))

    kept = len(edges)
    print(f"[stats] rows read={rows_total:,}  edges kept={kept:,}  "
          f"nodes={len(nodes):,}  self-loops dropped={self_loops:,}")
    print(f"[stats] interaction types: {type_counts}")

    write_provenance(config.PROCESSED_DIR / f"{LAYER}.meta.json", {
        "layer": LAYER,
        "source_url": URL,
        "biogrid_version": version,
        "filters": {"taxon": config.HUMAN_TAXON, "drop_self_loops": drop_self_loops},
        "rows_read": rows_total,
        "edges": kept,
        "nodes": len(nodes),
        "interaction_types": type_counts,
        "output": out_path.name,
        "detail": detail_path.name,
    })
    print(f"[ok] wrote {out_path}")
    print(f"[ok] wrote {detail_path}")
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(description="Build the BioGRID layer.")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    ap.add_argument("--keep-self-loops", action="store_true",
                    help="keep gene self-interactions (default: drop)")
    args = ap.parse_args()
    build(force=args.force, drop_self_loops=not args.keep_self_loops)


if __name__ == "__main__":
    main()
