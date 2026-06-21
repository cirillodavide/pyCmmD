"""Layer 3 -- Monarch (MONDO) shared-disease associations.

Connects two human genes when both are associated with the same disease (MONDO),
using the Monarch Initiative knowledge graph gene->disease edges.

Monarch's current KG sources these from OMIM / Orphanet / ClinGen only
(no GWAS). Associations split into:
  * causal      (biolink:CausalGeneToDiseaseAssociation)
  * correlated  (biolink:CorrelatedGeneToDiseaseAssociation)

Edge semantics (a gene pair derives from two gene->disease links):
  * interaction = "causal" | "correlated" | "mixed"  (mixed if the two differ)
  * annotation  = the disease name connecting the two genes
Nodes are Official Gene Symbols; the stable HGNC id is kept as source_id/target_id.
A companion detail file records, per edge, the MONDO id and each gene's
association category, predicate(s) and knowledge source(s).

Run:  python -m pycmmd.layers.monarch [--force]
"""

from __future__ import annotations

import argparse
import csv
import io
import tarfile
from collections import defaultdict
from itertools import combinations
from pathlib import Path

from .. import config
from ..utils import download, write_provenance

LAYER = "Monarch"
URL = "https://data.monarchinitiative.org/monarch-kg/latest/monarch-kg.tar.gz"

CATEGORY_MAP = {
    "biolink:CausalGeneToDiseaseAssociation": "causal",
    "biolink:CorrelatedGeneToDiseaseAssociation": "correlated",
}


def _stream(tf: tarfile.TarFile, member: str):
    return io.TextIOWrapper(tf.extractfile(member), encoding="utf-8", newline="")


def build(force: bool = False) -> Path:
    config.ensure_dirs()
    tar_path = download(URL, config.RAW_DIR / "monarch-kg.tar.gz", force=force, timeout=600)

    tf = tarfile.open(tar_path)

    # --- pass 1: collect gene<->disease associations ------------------------
    # disease -> gene_hgnc -> {"cats": set, "preds": set, "srcs": set}
    dz: dict[str, dict] = defaultdict(lambda: defaultdict(
        lambda: {"cats": set(), "preds": set(), "srcs": set()}))
    print("[parse] scanning edges for gene-disease associations ...")
    f = _stream(tf, "monarch-kg_edges.tsv")
    H = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    ip, ic, iks, isub, iobj = H["predicate"], H["category"], H["primary_knowledge_source"], H["subject"], H["object"]
    n_assoc = 0
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= iobj:
            continue
        cat = CATEGORY_MAP.get(p[ic])
        if not cat:
            continue
        sub, obj = p[isub], p[iobj]
        if not (sub.startswith("HGNC:") and obj.startswith("MONDO:")):
            continue
        g = dz[obj][sub]
        g["cats"].add(cat)
        g["preds"].add(p[ip].split(":")[-1])
        g["srcs"].add(p[iks].replace("infores:", ""))
        n_assoc += 1
    print(f"[parse] {n_assoc:,} gene-disease associations across {len(dz):,} diseases")

    needed_genes = {g for genes in dz.values() for g in genes}
    needed_dis = set(dz)

    # --- pass 2: id -> label maps from nodes --------------------------------
    print("[parse] resolving HGNC symbols and MONDO names from nodes ...")
    f = _stream(tf, "monarch-kg_nodes.tsv")
    N = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    nid, nname, nsym = N["id"], N["name"], N["symbol"]
    sym: dict[str, str] = {}
    dname: dict[str, str] = {}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= nname:
            continue
        i = p[nid]
        if i in needed_genes:
            sym[i] = (p[nsym] if len(p) > nsym and p[nsym] else p[nname]) or i
        elif i in needed_dis:
            dname[i] = p[nname] or i

    def gene_cat(info):  # one gene's category for a disease
        return "causal" if "causal" in info["cats"] else "correlated"

    # --- build edges --------------------------------------------------------
    out_path = config.PROCESSED_DIR / f"{LAYER}.edges.tsv"
    detail_path = config.PROCESSED_DIR / f"{LAYER}.edges.detail.tsv"
    seen: set = set()
    nodes: set = set()
    interaction_counts: dict[str, int] = {}

    with open(out_path, "w", newline="", encoding="utf-8") as out, \
         open(detail_path, "w", newline="", encoding="utf-8") as det:
        w = csv.writer(out, delimiter="\t")
        dw = csv.writer(det, delimiter="\t")
        w.writerow(config.EDGE_COLUMNS)
        dw.writerow(["source", "target", "mondo_id", "disease_name",
                     "source_category", "target_category",
                     "source_predicates", "target_predicates",
                     "source_sources", "target_sources"])
        for disease, genes in dz.items():
            if len(genes) < 2:
                continue
            label = dname.get(disease, disease)
            items = sorted(genes.items())  # by HGNC id, deterministic
            for (ha, ia), (hb, ib) in combinations(items, 2):
                sa, sb = sym.get(ha, ha), sym.get(hb, hb)
                if sa == sb:
                    continue
                if sb < sa:  # order by symbol; keep ids aligned
                    sa, sb, ha, hb, ia, ib = sb, sa, hb, ha, ib, ia
                key = (sa, sb, disease)
                if key in seen:
                    continue
                seen.add(key)
                ca, cb = gene_cat(ia), gene_cat(ib)
                interaction = ca if ca == cb else "mixed"
                w.writerow([sa, sb, LAYER, interaction, label, ha, hb])
                dw.writerow([sa, sb, disease, label, ca, cb,
                             ";".join(sorted(ia["preds"])), ";".join(sorted(ib["preds"])),
                             ";".join(sorted(ia["srcs"])), ";".join(sorted(ib["srcs"]))])
                nodes.update((sa, sb))
                interaction_counts[interaction] = interaction_counts.get(interaction, 0) + 1

    shared = sum(1 for genes in dz.values() if len(genes) >= 2)
    print(f"[stats] edges={len(seen):,}  nodes={len(nodes):,}  "
          f"diseases linking >=2 genes={shared:,}")
    print(f"[stats] interaction types: {interaction_counts}")

    write_provenance(config.PROCESSED_DIR / f"{LAYER}.meta.json", {
        "layer": LAYER,
        "source_url": URL,
        "gene_disease_associations": n_assoc,
        "diseases_total": len(dz),
        "diseases_linking_genes": shared,
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
    ap = argparse.ArgumentParser(description="Build the Monarch shared-disease layer.")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    args = ap.parse_args()
    build(force=args.force)


if __name__ == "__main__":
    main()
