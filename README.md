# pyCmmD

**An end-to-end Python pipeline to build a human gene *multilayer* network from public databases and analyze it with CmmD (Community multilayer Detection across modularity resolutions).**

pyCmmD downloads gene-association data from five public sources, turns each into a standardized gene–gene network layer, harmonizes all gene identifiers to a single HGNC-based index, and runs a multilayer community-detection analysis that tracks how persistently genes co-cluster as the modularity resolution is varied. The end products are **two matrices**:

1. a **community-membership matrix** — each gene's community at every resolution, and
2. a **pairwise Hamming-distance matrix** — how persistently each pair of genes stays in the same community across resolutions (low distance = robust co-membership).

This is a modern, pure-Python reimplementation of the workflow in
[cirillodavide/gene_multilayer_network](https://github.com/cirillodavide/gene_multilayer_network),
with refreshed data sources, identifier harmonization, and a parallelized CmmD engine.

---

## Scope

- **Organism:** human only (NCBI taxon 9606).
- **Nodes:** human genes, labeled by HGNC-approved gene symbol (stable `hgnc_id` retained alongside).
- **Edges:** gene–gene relationships, one *layer* per biological relationship type.
- **Analysis:** multilayer community detection (CmmD) across a sweep of modularity resolutions.
- **Not in scope:** non-human species, transcript/protein-isoform–level modeling, or directed/temporal networks.

---

## The five layers

Every layer is written in one common schema so they stack cleanly:

```
source · target · layer · interaction · annotation · source_id · target_id
```
where `source`/`target` are gene symbols, `interaction` is the controlled relationship type, `annotation` is the finer evidence (the connecting entity or assay), and `source_id`/`target_id` are the source database's native stable identifier.

| Layer | Relationship | Source (live, latest) | `interaction` types | Native id |
|-------|--------------|-----------------------|---------------------|-----------|
| **BioGRID** | molecular interaction (A binds/affects B) | BioGRID tab3 release | `physical`, `genetic` | Entrez |
| **ChEMBL** | drug co-targeting (same drug hits both genes) | ChEMBL mechanism-of-action API | `co-targeted` | UniProt |
| **Monarch** | shared disease (both genes ↔ same MONDO disease) | Monarch KG | `causal`, `correlated`, `mixed` | HGNC |
| **Reactome** | shared pathway (co-membership) | Reactome `NCBI2Reactome` | `TAS`, `IEA`, `mixed` | Entrez |
| **Recon3D** | shared metabolite (metabolic flow) | BiGG Recon3D model (JSON) | `producer-consumer`, `co-producer`, `co-consumer` | Entrez |

**Per-edge evidence** is preserved in companion `*.edges.detail.tsv` files (e.g. BioGRID PubMed IDs/throughput, ChEMBL drug+phase+mechanism, Monarch MONDO id + per-gene predicate/source, Reactome pathway id + evidence, Recon3D metabolite id + role).

**Co-membership layers** (Monarch, Reactome, Recon3D) turn each group (disease / pathway / metabolite) into a clique among its genes; a few very large groups can dominate, so each offers a size cap (`max_pathway_size`, `max_metabolite_size`). Recon3D additionally prunes super-connected "currency" metabolites (ATP, H2O, NAD, …) via [`pycmmd/resources/metabolites_to_prune.txt`](pycmmd/resources/metabolites_to_prune.txt).

---

## Harmonization to a common gene index

Different databases identify genes differently (Entrez / UniProt / HGNC / symbols / synonyms). The harmonization step ([`pycmmd/harmonize.py`](pycmmd/harmonize.py)) maps every node to a single canonical **(HGNC approved symbol, `hgnc_id`)** using the [HGNC complete set](https://www.genenames.org/download/), resolving in priority order:

1. the layer's **native stable id** (Entrez / UniProt / HGNC id),
2. the **approved symbol**,
3. an **unambiguous** previous/alias symbol (ambiguous synonyms are discarded, never guessed).

- Canonical identifiers are **unique and 1:1** with `hgnc_id`.
- Unresolvable nodes (obsolete IDs, withdrawn symbols) are **kept under their original symbol, flagged** with an empty `hgnc_id` (lossless) and listed per layer in `*.unmapped.tsv`.
- A unified gene index is written to `harmonized/nodes.tsv`.

---

## The CmmD analysis

CmmD ([`pycmmd/cmmd.py`](pycmmd/cmmd.py)) detects communities **jointly across all layers** for a sweep of modularity resolutions, then measures how stable each gene pair's co-membership is.

- **Engine:** multiplex modularity optimization via [`leidenalg`](https://github.com/vtraag/leidenalg) (`RBConfigurationVertexPartition`). The original CmmD drives MolTi; this is a pure-Python equivalent with the same objective. To match MolTi's *sum of per-layer normalized modularities*, per-layer weights default to `1/(2·mₐ)` (`layer_weighting = "normalized"`).
- **Resolution sweep:** `resolution_start … resolution_end` in steps of `interval` (paper default `0 → 30` step `0.5`).
- **Membership matrix:** genes × resolutions; each cell is the community id, plus a concatenated `Pattern` column. Two genes with the same `Pattern` were never separated across the whole sweep.
- **Hamming distance:** for each gene pair, the fraction of resolutions in which they fall in *different* communities. Stored **compressed** as `uint8` disagreement counts (decode: `distance = counts / n_resolutions`) — ~4× smaller than float32 and exact.
- **Parallelism:** both the resolution sweep and the distance computation are parallelized across processes (`jobs`).

> **Note on `n_iterations`:** `-1` runs Leiden to convergence (slower, fully stable partitions — use for final results); `2` (leidenalg's default) is much faster with near-identical partitions — good for exploration.

---

## Installation

Requires **Python ≥ 3.11** (uses the stdlib `tomllib`).

```bash
git clone <your-repo-url> pyCmmD
cd pyCmmD
python -m pip install -r requirements.txt
```

Dependencies: `numpy`, `scipy`, `pandas`, `python-igraph`, `leidenalg` (see [`requirements.txt`](requirements.txt)). No R or external binaries needed.

---

## Quick start

Run the entire pipeline (download → 5 layers → harmonize → CmmD) into a fresh dated session folder:

```bash
python -m pycmmd.run_all
```

This reads [`config.toml`](config.toml) and writes everything to `data/runs/<timestamp>/`. Source downloads (~700 MB total) are cached in `data/raw/` and shared across runs. A first run re-downloads everything; later runs reuse the cache unless `force_download = true`.

---

## Configuration

All tunable parameters live in [`config.toml`](config.toml):

```toml
[run]
jobs = 6                 # parallel worker processes for CmmD
force_download = false   # re-download source files even if cached

[chembl]
min_phase = 1.0          # 1 = approved + clinical, 4 = approved only

[reactome]
max_pathway_size = 0     # 0 = no cap; e.g. 100 drops huge cliques

[recon3d]
max_metabolite_size = 0  # 0 = no cap

[cmmd]
resolution_start = 0.0
resolution_end   = 30.0
interval         = 0.5
seed             = 42
layer_weighting  = "normalized"   # or "equal"
layers           = []             # which layers to use; [] = all (see "Choosing layers")
n_iterations     = -1             # -1 = converge (slow); 2 = fast
nodelist         = ""             # restrict the distance matrix to these genes; "" = all
compute_distance = true
```

### Choosing which layers to analyze

CmmD uses all five layers by default, but you can run on any subset (the node
universe becomes the union of the selected layers). This affects CmmD only — all
layers are still built and harmonized, so you can try combinations without
re-downloading. Set it in `config.toml` (`layers = ["BioGRID", "ChEMBL", "Reactome"]`)
or on the CLI (`--layers BioGRID ChEMBL Reactome`). This is also the easy way to
drop layers with commercial-use restrictions (Recon3D; Monarch via OMIM) — see
**Data sources** for licensing.

Each run **copies the config it used** into its session folder, and records per-step timings in `manifest.json`, for full reproducibility.

---

## Output layout

```
data/
├── raw/                         # shared download cache (BioGRID, Monarch KG, Recon3D, HGNC, …)
└── runs/<timestamp>/
    ├── config.toml              # the exact config used (provenance)
    ├── manifest.json            # per-step timings + total
    ├── processed/               # one standardized layer per source
    │   ├── <Layer>.edges.tsv          # common schema
    │   ├── <Layer>.edges.detail.tsv   # full per-edge evidence
    │   └── <Layer>.meta.json          # provenance (resolved DB version, filters, counts)
    ├── harmonized/
    │   ├── <Layer>.edges.tsv          # nodes remapped to HGNC; source_id/target_id = hgnc_id
    │   ├── <Layer>.unmapped.tsv       # flagged unresolved identifiers
    │   ├── nodes.tsv                  # unified gene index (symbol, hgnc_id, harmonized, layers)
    │   └── harmonization_report.json
    └── cmmd/
        ├── membership_matrix.tsv      # genes × resolutions (+ Pattern)   ← matrix 1
        ├── membership.npy             # same, as int array
        ├── hamming_counts.npy         # compressed uint8 condensed        ← matrix 2
        ├── hamming_nodes.txt          # gene order for the distance matrix
        └── cmmd_report.json           # params, communities/resolution, decode info
```

---

## Working with the results

**Look up distances without loading the whole matrix** (memory-mapped accessor):

```python
from pycmmd.distance import CmmDDistance
d = CmmDDistance("data/runs/<timestamp>/cmmd")   # or default data/cmmd/
d.dist("BRCA1", "BARD1")     # one decoded distance
d.row("TP53")                # distances to all genes
d.nearest("TP53", k=10)      # most persistently co-clustered genes
```

**Materialize the full distance matrix** from the compressed counts:

```bash
python -m pycmmd.distance --form square        # → hamming_square.npy (~1.7 GB float32)
python -m pycmmd.distance --form condensed     # → condensed float32 (~870 MB)
python -m pycmmd.distance --form square --tsv  # also a labeled TSV (large/slow)
```
Decode rule: **`distance = hamming_counts / n_resolutions`** (the divisor is stored in `cmmd_report.json`).

---

## Running steps individually

For development/debugging, each stage is runnable on its own (writing to `data/<name>/` by default):

```bash
python -m pycmmd.layers.biogrid      [--force] [--keep-self-loops]
python -m pycmmd.layers.chembl       [--min-phase 4]
python -m pycmmd.layers.monarch      [--force]
python -m pycmmd.layers.reactome     [--force] [--max-pathway-size 100]
python -m pycmmd.layers.recon3d      [--force] [--max-metabolite-size 100]
python -m pycmmd.harmonize           [--force]
python -m pycmmd.cmmd                [--resolution-end 30 --interval 0.5 --jobs 6 --n-iterations 2 --layers BioGRID ChEMBL Reactome --no-distance]
```

---

## Reproducibility

- Each run is an isolated, timestamped session — **re-running never overwrites** previous results.
- Every layer records its **resolved database version** in `*.meta.json`; the raw downloads are the version snapshot (cached in `data/raw/`).
- CmmD is **deterministic** given the `seed`.
- The config used is archived inside each session.

> Because sources are pulled "latest" by default, results will evolve as the underlying databases update. Pin to a `data/raw/` snapshot (and keep `force_download = false`) to reproduce a specific run.

---

## Data sources

| Source | URL | Notes |
|--------|-----|-------|
| BioGRID | https://thebiogrid.org | molecular interactions (tab3) |
| ChEMBL | https://www.ebi.ac.uk/chembl | drug mechanism of action (REST API) |
| Monarch Initiative | https://monarchinitiative.org | gene–disease (knowledge graph) |
| Reactome | https://reactome.org | pathway annotations |
| BiGG / Recon3D | http://bigg.ucsd.edu | human metabolic model |
| HGNC | https://www.genenames.org | gene nomenclature / identifier harmonization |
| NCBI Gene | https://www.ncbi.nlm.nih.gov/gene | Entrez↔symbol mapping (`gene_info`) |

Please cite the individual databases and the original
[CmmD method](https://github.com/ikernunezca/CmmD) when using results from this pipeline.
All source data remain subject to the licenses of their respective providers.

---

## Project structure

```
pycmmd/
├── config.py           # paths, common schema, dated-session redirection
├── utils.py            # caching downloader, JSON fetch, provenance writer
├── mappings.py         # NCBI gene_info Entrez→Symbol helper
├── layers/             # one builder per source layer
│   ├── biogrid.py  chembl.py  monarch.py  reactome.py  recon3d.py
├── harmonize.py        # HGNC-based common gene index
├── cmmd.py             # multilayer community detection + Hamming (parallel)
├── distance.py         # lazy distance accessor + full-matrix export CLI
├── run_all.py          # end-to-end orchestrator (reads config.toml)
└── resources/
    └── metabolites_to_prune.txt
config.toml             # all pipeline parameters
requirements.txt
```

---

## License

The original CmmD method is GPL (≥2); the source databases carry their own licenses.
