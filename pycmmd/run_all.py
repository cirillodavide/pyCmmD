"""End-to-end pipeline orchestrator.

Reads config.toml, then runs the whole pipeline in one process into a fresh
dated session folder (data/runs/<timestamp>/): 5 layers -> harmonize -> CmmD.
Nothing from previous runs is overwritten. The raw download cache (data/raw/)
is shared across runs; each run records the DB versions it used in its meta.

Run:  python -m pycmmd.run_all [--config config.toml] [--session-name NAME]
"""

from __future__ import annotations

import argparse
import shutil
import time
import tomllib
from datetime import datetime
from pathlib import Path

from . import cmmd, config, harmonize
from .layers import biogrid, chembl, monarch, reactome, recon3d
from .utils import write_provenance

DEFAULT_CONFIG = config.ROOT / "config.toml"


def load_config(path: Path) -> dict:
    if path and Path(path).exists():
        with open(path, "rb") as fh:
            return tomllib.load(fh)
    print(f"[run_all] no config at {path}; using built-in defaults")
    return {}


def _cap(value):
    """Treat 0 / empty as 'no cap' (None)."""
    return value or None


def main() -> None:
    ap = argparse.ArgumentParser(description="Run the full pyCmmD pipeline into a dated session.")
    ap.add_argument("--config", default=str(DEFAULT_CONFIG), help="path to config.toml")
    ap.add_argument("--session-name", default=None, help="override the dated folder name")
    args = ap.parse_args()

    cfg = load_config(Path(args.config))
    run = cfg.get("run", {})
    force = bool(run.get("force_download", False))
    jobs = run.get("jobs", None)

    name = args.session_name or datetime.now().strftime("%Y-%m-%d_%H%M%S")
    session = config.RUNS_DIR / name
    config.use_session(session)
    print(f"[run_all] session -> {session}")
    if Path(args.config).exists():
        shutil.copy(args.config, session / "config.toml")

    bg = cfg.get("biogrid", {})
    ch = cfg.get("chembl", {})
    rx = cfg.get("reactome", {})
    rc = cfg.get("recon3d", {})
    cm = cfg.get("cmmd", {})

    steps = [
        ("BioGRID", lambda: biogrid.build(force=force,
                                          drop_self_loops=bg.get("drop_self_loops", True))),
        ("ChEMBL", lambda: chembl.build(min_phase=ch.get("min_phase", 1.0))),
        ("Monarch", lambda: monarch.build(force=force)),
        ("Reactome", lambda: reactome.build(force=force,
                                            max_pathway_size=_cap(rx.get("max_pathway_size", 0)))),
        ("Recon3D", lambda: recon3d.build(force=force,
                                          max_metabolite_size=_cap(rc.get("max_metabolite_size", 0)))),
        ("Harmonize", lambda: harmonize.build(force=force)),
        ("CmmD", lambda: cmmd.run_cmmd(
            resolution_start=cm.get("resolution_start", 0.0),
            resolution_end=cm.get("resolution_end", 30.0),
            interval=cm.get("interval", 0.5),
            seed=cm.get("seed", 42),
            layer_weighting=cm.get("layer_weighting", "normalized"),
            n_iterations=cm.get("n_iterations", -1),
            nodelist=(cm.get("nodelist", "") or None),
            layers=(cm.get("layers") or None),
            compute_distance=cm.get("compute_distance", True),
            jobs=jobs)),
    ]

    t0 = time.time()
    timings = {}
    for label, fn in steps:
        print(f"\n========== {label} ==========")
        s = time.time()
        fn()
        timings[label] = round(time.time() - s, 1)
        print(f"[run_all] {label} done in {timings[label]:.1f}s")

    total = round(time.time() - t0, 1)
    write_provenance(session / "manifest.json", {
        "session": name,
        "config": cfg,
        "step_seconds": timings,
        "total_seconds": total,
        "outputs": {"processed": "processed/", "harmonized": "harmonized/", "cmmd": "cmmd/"},
    })
    (config.RUNS_DIR / "latest.txt").write_text(str(session), encoding="utf-8")

    print(f"\n[run_all] COMPLETE in {total/60:.1f} min -> {session}")
    print("[run_all] per-step seconds:", timings)


if __name__ == "__main__":
    main()
