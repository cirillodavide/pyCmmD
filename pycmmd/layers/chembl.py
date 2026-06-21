"""Layer 2 -- ChEMBL drug co-targeting.

Connects two human genes when the same drug targets both (or their products),
using ChEMBL's curated *mechanism of action* table as the drug->target source.

Scope: drugs at max_phase >= 1 (approved + clinical); the actual phase is
recorded per edge.

Edge semantics:
  * interaction = "co-targeted"  (constant for this layer)
  * annotation  = the drug name connecting the two genes
Nodes are Official Gene Symbols; UniProt accession is kept as source_id/target_id.
A companion detail file preserves drug ChEMBL id, phase, action type and the
mechanism-of-action text per edge.

Run:  python -m pycmmd.layers.chembl [--min-phase N]
"""

from __future__ import annotations

import argparse
import csv
from itertools import combinations
from urllib.parse import quote

from .. import config
from ..utils import get_json, write_provenance

LAYER = "ChEMBL"
API = "https://www.ebi.ac.uk/chembl/api/data"
PAGE = 1000      # mechanism page size (ChEMBL caps pages at 1000)
BATCH = 40       # ids per __in query


def _paged(endpoint: str, root_key: str):
    """Yield every record from a paginated ChEMBL endpoint."""
    url = f"{API}/{endpoint}.json?limit={PAGE}"
    while url:
        data = get_json(url, retries=8)
        yield from data[root_key]
        nxt = data["page_meta"]["next"]
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def fetch_mechanisms(min_phase: float):
    """Return mechanism rows with max_phase >= min_phase."""
    rows = []
    for m in _paged("mechanism", "mechanisms"):
        phase = m.get("max_phase")
        if phase is None or float(phase) < min_phase:
            continue
        if not m.get("target_chembl_id") or not m.get("molecule_chembl_id"):
            continue
        rows.append({
            "molecule": m["molecule_chembl_id"],
            "target": m["target_chembl_id"],
            "action_type": m.get("action_type") or "UNKNOWN",
            "max_phase": phase,
            "moa": m.get("mechanism_of_action") or "",
        })
    return rows


def fetch_target_genes(target_ids):
    """Map human target_chembl_id -> list of (gene_symbol, uniprot)."""
    out = {}
    ids = sorted(target_ids)
    for batch in _chunks(ids, BATCH):
        q = quote(",".join(batch))
        data = get_json(f"{API}/target.json?target_chembl_id__in={q}&limit=1000", retries=8)
        for t in data["targets"]:
            if t.get("organism") != "Homo sapiens":
                continue
            genes = []
            for comp in t.get("target_components") or []:
                sym = next((s["component_synonym"]
                            for s in comp.get("target_component_synonyms") or []
                            if s.get("syn_type") == "GENE_SYMBOL"), None)
                if sym:
                    genes.append((sym.strip(), (comp.get("accession") or "").strip()))
            if genes:
                out[t["target_chembl_id"]] = genes
    return out


def fetch_drug_names(mol_ids):
    """Map molecule_chembl_id -> preferred drug name (fallback to the id)."""
    out = {}
    ids = sorted(mol_ids)
    for batch in _chunks(ids, BATCH):
        q = quote(",".join(batch))
        data = get_json(f"{API}/molecule.json?molecule_chembl_id__in={q}&limit=1000", retries=8)
        for mol in data["molecules"]:
            cid = mol["molecule_chembl_id"]
            out[cid] = (mol.get("pref_name") or cid)
    return out


def build(min_phase: float = 1.0) -> str:
    config.ensure_dirs()

    print(f"[fetch] mechanism table (max_phase >= {min_phase}) ...")
    mechs = fetch_mechanisms(min_phase)
    print(f"[fetch] {len(mechs):,} mechanism rows kept")

    target_ids = {m["target"] for m in mechs}
    print(f"[fetch] resolving {len(target_ids):,} targets -> human genes ...")
    target_genes = fetch_target_genes(target_ids)
    print(f"[fetch] {len(target_genes):,} human targets resolved")

    # drug -> {gene_symbol: {"uniprot","action","moa"}} (human targets only)
    drug_targets: dict[str, dict] = {}
    phase_of: dict[str, float] = {}
    for m in mechs:
        genes = target_genes.get(m["target"])
        if not genes:
            continue
        d = drug_targets.setdefault(m["molecule"], {})
        phase_of[m["molecule"]] = max(phase_of.get(m["molecule"], 0), float(m["max_phase"]))
        for sym, uni in genes:
            g = d.setdefault(sym, {"uniprot": uni, "actions": set(), "moas": set()})
            g["actions"].add(m["action_type"])
            if m["moa"]:
                g["moas"].add(m["moa"])

    # only drugs hitting >= 2 distinct genes yield co-targeting edges
    edge_drugs = {d: g for d, g in drug_targets.items() if len(g) >= 2}
    print(f"[fetch] resolving names for {len(edge_drugs):,} multi-target drugs ...")
    names = fetch_drug_names(edge_drugs.keys())

    out_path = config.PROCESSED_DIR / f"{LAYER}.edges.tsv"
    detail_path = config.PROCESSED_DIR / f"{LAYER}.edges.detail.tsv"
    seen: set = set()
    nodes: set = set()
    phase_hist: dict[str, int] = {}

    with open(out_path, "w", newline="", encoding="utf-8") as out, \
         open(detail_path, "w", newline="", encoding="utf-8") as det:
        w = csv.writer(out, delimiter="\t")
        dw = csv.writer(det, delimiter="\t")
        w.writerow(config.EDGE_COLUMNS)
        dw.writerow(["source", "target", "drug_chembl_id", "drug_name",
                     "max_phase", "action_type", "mechanism_of_action"])
        for drug, genes in edge_drugs.items():
            name = names.get(drug, drug)
            phase = phase_of[drug]
            pkey = f"phase_{phase:g}"
            for (ga, ia), (gb, ib) in combinations(sorted(genes.items()), 2):
                sa, sb = ga, gb
                ua, ub = ia["uniprot"], ib["uniprot"]
                key = (sa, sb, drug)
                if key in seen:
                    continue
                seen.add(key)
                w.writerow([sa, sb, LAYER, "co-targeted", name, ua, ub])
                actions = ";".join(sorted(ia["actions"] | ib["actions"]))
                moas = " | ".join(sorted(ia["moas"] | ib["moas"]))
                dw.writerow([sa, sb, drug, name, f"{phase:g}", actions, moas])
                nodes.update((sa, sb))
                phase_hist[pkey] = phase_hist.get(pkey, 0) + 1

    print(f"[stats] edges={len(seen):,}  nodes={len(nodes):,}  "
          f"co-targeting drugs={len(edge_drugs):,}")
    print(f"[stats] edges by drug phase: {phase_hist}")

    write_provenance(config.PROCESSED_DIR / f"{LAYER}.meta.json", {
        "layer": LAYER,
        "source": f"{API}/mechanism",
        "min_phase": min_phase,
        "mechanism_rows": len(mechs),
        "human_targets": len(target_genes),
        "co_targeting_drugs": len(edge_drugs),
        "edges": len(seen),
        "nodes": len(nodes),
        "edges_by_phase": phase_hist,
        "output": out_path.name,
        "detail": detail_path.name,
    })
    print(f"[ok] wrote {out_path}")
    print(f"[ok] wrote {detail_path}")
    return str(out_path)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build the ChEMBL drug co-targeting layer.")
    ap.add_argument("--min-phase", type=float, default=1.0,
                    help="minimum drug max_phase (default 1 = approved+clinical; 4 = approved only)")
    args = ap.parse_args()
    build(min_phase=args.min_phase)


if __name__ == "__main__":
    main()
