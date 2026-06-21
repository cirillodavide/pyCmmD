"""Layer 5 -- Recon3D metabolic shared-metabolite associations.

Connects two human genes when their products participate in reactions sharing
the same metabolite, using the Recon3D human metabolic model (BiGG, JSON).
Super-connected currency metabolites (ATP, H2O, NAD, ...) are pruned.

Edge semantics (metabolite-mediated, like the original pipeline):
  * interaction = "producer-consumer" | "co-producer" | "co-consumer"
                  (from the metabolite's role -- reactant/product -- in each
                  gene's reactions; producer-consumer means metabolic flow)
  * annotation  = the shared metabolite name connecting the two genes
Nodes are Official Gene Symbols (carried by the model, Entrez->Symbol fallback
via NCBI gene_info); the stable Entrez id is kept as source_id/target_id.
A companion detail file records the metabolite id and each gene's role.

Run:  python -m pycmmd.layers.recon3d [--force] [--max-metabolite-size N]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from itertools import combinations
from pathlib import Path

from .. import config
from ..mappings import entrez_to_symbol
from ..utils import download, write_provenance

LAYER = "Recon3D"
URL = "http://bigg.ucsd.edu/static/models/Recon3D.json"
PRUNE_FILE = Path(__file__).resolve().parent.parent / "resources" / "metabolites_to_prune.txt"


def load_prune_set() -> set[str]:
    out = set()
    for line in PRUNE_FILE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.add(line)
    return out


def build(force: bool = False, max_metabolite_size: int | None = None) -> Path:
    config.ensure_dirs()
    src = download(URL, config.RAW_DIR / "Recon3D.json", force=force)
    model = json.loads(src.read_text(encoding="utf-8"))
    e2s = entrez_to_symbol(force=force)
    prune = load_prune_set()

    comps = set(model["compartments"])

    def strip_comp(mid: str) -> str:
        head, _, tail = mid.rpartition("_")
        return head if (head and tail in comps) else mid

    # gene_id -> (entrez, symbol); model carries the symbol as the gene name
    gene_info: dict[str, tuple[str, str]] = {}
    for g in model["genes"]:
        gid = g["id"]
        if gid == "0":
            continue
        entrez = gid.split("_")[0]
        symbol = (g.get("name") or "").strip() or e2s.get(entrez) or entrez
        gene_info[gid] = (entrez, symbol)
    gene_ids = set(gene_info)

    # metabolite base -> name (first seen)
    met_name: dict[str, str] = {}
    for x in model["metabolites"]:
        b = strip_comp(x["id"])
        met_name.setdefault(b, (x.get("name") or b).strip())

    token_re = re.compile(r"[()\s]+")

    def genes_of(rule: str):
        if not rule:
            return []
        return [t for t in token_re.split(rule) if t in gene_ids]

    # metabolite base -> gene_id -> set of roles ("reactant"/"product")
    met_genes: dict[str, dict] = defaultdict(lambda: defaultdict(set))
    for r in model["reactions"]:
        gs = genes_of(r.get("gene_reaction_rule", ""))
        if not gs:
            continue
        for mid, stoich in r["metabolites"].items():
            b = strip_comp(mid)
            if b in prune:
                continue
            role = "reactant" if stoich < 0 else "product"
            for gid in gs:
                met_genes[b][gid].add(role)

    def pair_interaction(r1: set, r2: set) -> str:
        flow = ("product" in r1 and "reactant" in r2) or ("reactant" in r1 and "product" in r2)
        if flow:
            return "producer-consumer"
        if r1 == {"product"} and r2 == {"product"}:
            return "co-producer"
        if r1 == {"reactant"} and r2 == {"reactant"}:
            return "co-consumer"
        return "producer-consumer"

    out_path = config.PROCESSED_DIR / f"{LAYER}.edges.tsv"
    detail_path = config.PROCESSED_DIR / f"{LAYER}.edges.detail.tsv"
    seen: set = set()
    nodes: set = set()
    interaction_counts: dict[str, int] = {}
    skipped_large = 0

    with open(out_path, "w", newline="", encoding="utf-8") as out, \
         open(detail_path, "w", newline="", encoding="utf-8") as det:
        w = csv.writer(out, delimiter="\t")
        dw = csv.writer(det, delimiter="\t")
        w.writerow(config.EDGE_COLUMNS)
        dw.writerow(["source", "target", "metabolite_id", "metabolite_name",
                     "source_role", "target_role"])
        for met, genes in met_genes.items():
            if len(genes) < 2:
                continue
            if max_metabolite_size and len(genes) > max_metabolite_size:
                skipped_large += 1
                continue
            name = met_name.get(met, met)
            # resolve to (symbol, entrez, roles), order by symbol
            members = sorted((gene_info[g][1], gene_info[g][0], r) for g, r in genes.items())
            for (sa, ea, ra), (sb, eb, rb) in combinations(members, 2):
                if sa == sb:
                    continue
                key = (sa, sb, met)
                if key in seen:
                    continue
                seen.add(key)
                interaction = pair_interaction(ra, rb)
                w.writerow([sa, sb, LAYER, interaction, name, ea, eb])
                dw.writerow([sa, sb, met, name, "+".join(sorted(ra)), "+".join(sorted(rb))])
                nodes.update((sa, sb))
                interaction_counts[interaction] = interaction_counts.get(interaction, 0) + 1

    print(f"[stats] reactions={len(model['reactions']):,}  genes={len(gene_ids):,}  "
          f"metabolites(pruned-applied)={len(met_genes):,}")
    print(f"[stats] edges={len(seen):,}  nodes={len(nodes):,}  "
          f"pruned currency metabolites={len(prune)}  large metabolites skipped={skipped_large}")
    print(f"[stats] interaction types: {interaction_counts}")

    write_provenance(config.PROCESSED_DIR / f"{LAYER}.meta.json", {
        "layer": LAYER,
        "source_url": URL,
        "model_id": model.get("id"),
        "model_version": model.get("version"),
        "reactions": len(model["reactions"]),
        "genes": len(gene_ids),
        "pruned_metabolites": sorted(prune),
        "max_metabolite_size": max_metabolite_size,
        "large_metabolites_skipped": skipped_large,
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
    ap = argparse.ArgumentParser(description="Build the Recon3D metabolic layer.")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    ap.add_argument("--max-metabolite-size", type=int, default=None,
                    help="exclude metabolites shared by more than N genes (damps hubs)")
    args = ap.parse_args()
    build(force=args.force, max_metabolite_size=args.max_metabolite_size)


if __name__ == "__main__":
    main()
