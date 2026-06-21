"""Layer 4 -- Reactome shared-pathway associations.

Connects two human genes when both belong to the same Reactome pathway
(lowest-level pathway annotations from NCBI2Reactome.txt).

Edge semantics (a gene pair derives from two gene->pathway annotations):
  * interaction = "TAS" | "IEA" | "mixed"   (Reactome evidence code; mixed if
                  the two genes' annotations differ)
  * annotation  = the pathway name connecting the two genes
Nodes are Official Gene Symbols (Entrez->Symbol via NCBI gene_info); the stable
Entrez id is kept as source_id/target_id. A companion detail file records the
R-HSA pathway id and each gene's evidence code.

Run:  python -m pycmmd.layers.reactome [--force] [--max-pathway-size N]
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from itertools import combinations
from pathlib import Path

from .. import config
from ..mappings import entrez_to_symbol
from ..utils import download, write_provenance

LAYER = "Reactome"
URL = "https://reactome.org/download/current/NCBI2Reactome.txt"

# NCBI2Reactome.txt columns (0-indexed):
# 0=EntrezGeneID 1=PathwayStableID 2=URL 3=PathwayName 4=EvidenceCode 5=Species
C_GENE, C_PATH, C_NAME, C_EVID, C_SPECIES = 0, 1, 3, 4, 5


def build(force: bool = False, max_pathway_size: int | None = None) -> Path:
    config.ensure_dirs()
    src = download(URL, config.RAW_DIR / "NCBI2Reactome.txt", force=force)
    e2s = entrez_to_symbol(force=force)

    # pathway -> {entrez: evidence_code}, plus pathway names
    pw: dict[str, dict] = defaultdict(dict)
    pw_name: dict[str, str] = {}
    rows = 0
    with open(src, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= C_SPECIES or p[C_SPECIES] != "Homo sapiens":
                continue
            rows += 1
            pw[p[C_PATH]][p[C_GENE]] = p[C_EVID]
            pw_name[p[C_PATH]] = p[C_NAME].strip()

    out_path = config.PROCESSED_DIR / f"{LAYER}.edges.tsv"
    detail_path = config.PROCESSED_DIR / f"{LAYER}.edges.detail.tsv"
    seen: set = set()
    nodes: set = set()
    interaction_counts: dict[str, int] = {}
    skipped_large = 0
    unmapped: set = set()

    with open(out_path, "w", newline="", encoding="utf-8") as out, \
         open(detail_path, "w", newline="", encoding="utf-8") as det:
        w = csv.writer(out, delimiter="\t")
        dw = csv.writer(det, delimiter="\t")
        w.writerow(config.EDGE_COLUMNS)
        dw.writerow(["source", "target", "pathway_id", "pathway_name",
                     "source_evidence", "target_evidence"])
        for pid, genes in pw.items():
            # resolve to symbols; keep (symbol, entrez, evidence)
            members = []
            for ent, ev in genes.items():
                s = e2s.get(ent)
                if s:
                    members.append((s, ent, ev))
                else:
                    unmapped.add(ent)
            if len(members) < 2:
                continue
            if max_pathway_size and len(members) > max_pathway_size:
                skipped_large += 1
                continue
            name = pw_name.get(pid, pid)
            for (sa, ea, va), (sb, eb, vb) in combinations(sorted(members), 2):
                if sa == sb:
                    continue
                key = (sa, sb, pid)
                if key in seen:
                    continue
                seen.add(key)
                interaction = va if va == vb else "mixed"
                w.writerow([sa, sb, LAYER, interaction, name, ea, eb])
                dw.writerow([sa, sb, pid, name, va, vb])
                nodes.update((sa, sb))
                interaction_counts[interaction] = interaction_counts.get(interaction, 0) + 1

    print(f"[stats] human annotations={rows:,}  pathways={len(pw):,}")
    print(f"[stats] edges={len(seen):,}  nodes={len(nodes):,}  "
          f"unmapped Entrez={len(unmapped):,}  large pathways skipped={skipped_large}")
    print(f"[stats] interaction types: {interaction_counts}")

    write_provenance(config.PROCESSED_DIR / f"{LAYER}.meta.json", {
        "layer": LAYER,
        "source_url": URL,
        "mapping": "NCBI Homo_sapiens.gene_info (Entrez->Symbol)",
        "human_annotations": rows,
        "pathways": len(pw),
        "max_pathway_size": max_pathway_size,
        "large_pathways_skipped": skipped_large,
        "unmapped_entrez": len(unmapped),
        "edges": len(seen),
        "nodes": len(nodes),
        "interaction_types": interaction_counts,
        "output": out_path.name,
        "detail": detail_path.name,
    })
    print(f"[ok] wrote {out_path}")
    print(f"[ok] wrote {detail_path}")
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(description="Build the Reactome shared-pathway layer.")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    ap.add_argument("--max-pathway-size", type=int, default=None,
                    help="exclude pathways with more than N genes (damps huge cliques)")
    args = ap.parse_args()
    build(force=args.force, max_pathway_size=args.max_pathway_size)


if __name__ == "__main__":
    main()
