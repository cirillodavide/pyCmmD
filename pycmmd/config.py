"""Central configuration: paths and the common edge schema shared by all layers."""

from __future__ import annotations

from pathlib import Path

# Project layout -------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"          # downloaded archives / source files (shared cache)
RUNS_DIR = DATA_DIR / "runs"        # dated session folders (one per full run_all)

# Derived-output dirs. Default to data/<name> for standalone module runs; a full
# pipeline run redirects them into a dated session folder via use_session().
PROCESSED_DIR = DATA_DIR / "processed"    # standardized per-layer edge lists
HARMONIZED_DIR = DATA_DIR / "harmonized"  # per-layer edges remapped to HGNC index
CMMD_DIR = DATA_DIR / "cmmd"        # CmmD membership + Hamming distance outputs

# Human NCBI Taxonomy identifier.
HUMAN_TAXON = "9606"

# Common edge schema shared by every layer -----------------------------------
# Each layer writes a TSV with exactly these columns so they can be stacked
# into one multilayer edge list and harmonized to a common index later.
#
#   source, target : node identifiers (Official Gene Symbols for now)
#   layer          : the layer name (e.g. "BioGRID")
#   interaction    : the controlled interaction *type* (e.g. physical/genetic)
#   annotation     : finer-grained evidence/detail (e.g. the assay used)
#   source_id,
#   target_id      : the native, stable database identifier when available
#                    (e.g. Entrez Gene ID) -- kept to anchor harmonization.
EDGE_COLUMNS = [
    "source",
    "target",
    "layer",
    "interaction",
    "annotation",
    "source_id",
    "target_id",
]


def ensure_dirs() -> None:
    """Create the data directories if they do not exist."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    HARMONIZED_DIR.mkdir(parents=True, exist_ok=True)
    CMMD_DIR.mkdir(parents=True, exist_ok=True)


def use_session(session_dir) -> Path:
    """Redirect derived outputs (processed/harmonized/cmmd) into a session dir.

    The raw download cache stays shared at data/raw/. Called by run_all so a full
    pipeline run lands in one dated folder without overwriting previous runs.
    """
    global PROCESSED_DIR, HARMONIZED_DIR, CMMD_DIR
    session_dir = Path(session_dir)
    PROCESSED_DIR = session_dir / "processed"
    HARMONIZED_DIR = session_dir / "harmonized"
    CMMD_DIR = session_dir / "cmmd"
    ensure_dirs()
    return session_dir
