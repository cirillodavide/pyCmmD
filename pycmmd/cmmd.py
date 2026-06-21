"""CmmD -- Community multilayer Detection across modularity resolutions.

Pure-Python reimplementation of the CmmD algorithm (Nunez-Carpintero et al.),
which originally drives MolTi. For each resolution gamma in a sweep, it runs
multiplex modularity community detection jointly over all layers (leidenalg,
RBConfigurationVertexPartition), then tracks each gene's community across
resolutions and measures persistence of co-membership via Hamming distance.

Faithfulness note: MolTi maximizes the SUM of each layer's (normalized)
modularity. leidenalg's RBConfigurationVertexPartition quality is unnormalized,
so per-layer weights default to 1/(2*m_layer) ("normalized"), matching MolTi's
objective. Optimizer is Leiden (vs MolTi's Louvain): same objective.

Performance:
  * the resolution sweep is parallelized across processes (independent gammas)
  * the Hamming distance is parallelized across processes, writing into one
    shared-memory buffer
Compression: every Hamming value is k/R for an integer disagreement count
k in 0..R, so the distance is stored as uint8 counts (4x smaller than float32,
exact, memory-mappable). Decode with distance = counts / n_resolutions; the
``pycmmd.distance.CmmDDistance`` accessor does this lazily.

Outputs (data/cmmd/):
  * membership_matrix.tsv / membership.npy  -- genes x resolutions community ids
  * hamming_counts.npy (uint8 condensed) + hamming_nodes.txt
  * cmmd_report.json  -- params, per-resolution communities, decode info

Run:  python -m pycmmd.cmmd [--resolution-start 0 --resolution-end 30 --interval 0.5]
                            [--jobs 6] [--nodelist genes.txt] [--seed 42] [--no-distance]
"""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import shared_memory
from pathlib import Path

import igraph as ig
import leidenalg as la
import numpy as np

from . import config
from .harmonize import LAYER_NATIVE
from .utils import write_provenance


# --------------------------------------------------------------------------- #
# Multiplex loading
# --------------------------------------------------------------------------- #
def load_multiplex(layer_dir: Path):
    """Load harmonized layers into igraphs over a shared vertex set.

    Returns (layers, graphs, node_names). Each layer is collapsed to simple
    undirected edges; all graphs share vertex indices 0..n-1.
    """
    layer_edges: dict[str, set] = {}
    nodes: set[str] = set()
    for layer in LAYER_NATIVE:
        path = layer_dir / f"{layer}.edges.tsv"
        if not path.exists():
            continue
        pairs = set()
        with open(path, encoding="utf-8") as fh:
            r = csv.reader(fh, delimiter="\t")
            next(r)
            for row in r:
                a, b = row[0], row[1]
                pairs.add((a, b) if a < b else (b, a))
                nodes.add(a)
                nodes.add(b)
        layer_edges[layer] = pairs

    names = sorted(nodes)
    idx = {n: i for i, n in enumerate(names)}
    layers, graphs = [], []
    for layer, pairs in layer_edges.items():
        g = ig.Graph(n=len(names),
                     edges=[(idx[a], idx[b]) for a, b in pairs], directed=False)
        layers.append(layer)
        graphs.append(g)
    return layers, graphs, names


def _layer_weights(graphs, scheme: str):
    if scheme == "equal":
        return [1.0] * len(graphs)
    return [1.0 / (2.0 * g.ecount()) if g.ecount() else 0.0 for g in graphs]


# --------------------------------------------------------------------------- #
# Parallel resolution sweep
# --------------------------------------------------------------------------- #
_SW: dict = {}


def _init_sweep(harmonized_dir: str, weighting: str, seed: int, n_iterations: int):
    # harmonized_dir passed explicitly so spawned workers don't depend on the
    # parent's (possibly session-redirected) config globals.
    _, graphs, names = load_multiplex(Path(harmonized_dir))
    _SW["graphs"] = graphs
    _SW["weights"] = _layer_weights(graphs, weighting)
    _SW["seed"] = seed
    _SW["n_iterations"] = n_iterations
    _SW["n"] = len(names)


def _sweep_one(gamma: float):
    graphs, weights, seed = _SW["graphs"], _SW["weights"], _SW["seed"]
    optimiser = la.Optimiser()
    optimiser.set_rng_seed(seed)
    parts = [la.RBConfigurationVertexPartition(g, resolution_parameter=gamma) for g in graphs]
    optimiser.optimise_partition_multiplex(parts, layer_weights=weights,
                                           n_iterations=_SW["n_iterations"])
    return np.asarray(parts[0].membership, dtype=np.int32)


def sweep(resolutions, harmonized_dir, weighting, seed, jobs, n_iterations):
    """Return (membership n x R int32, communities-per-resolution list)."""
    initargs = (str(harmonized_dir), weighting, seed, n_iterations)
    if jobs <= 1:
        _init_sweep(*initargs)
        cols = [_sweep_one(g) for g in resolutions]
    else:
        with ProcessPoolExecutor(max_workers=jobs, initializer=_init_sweep,
                                 initargs=initargs) as ex:
            cols = list(ex.map(_sweep_one, resolutions))
    membership = np.column_stack(cols).astype(np.int32)
    communities = [int(c.max()) + 1 for c in cols]
    return membership, communities


# --------------------------------------------------------------------------- #
# Parallel Hamming distance (disagreement counts, uint8/uint16)
# --------------------------------------------------------------------------- #
_HM: dict = {}


def _row_offset(i: int, n: int) -> int:
    """Start index of row i in the condensed upper triangle."""
    return n * (n - 1) // 2 - (n - i) * (n - i - 1) // 2


def _init_ham(M, shm_name, n, R, dtype_str, npairs):
    _HM.update(M=M, shm_name=shm_name, n=n, R=R, dtype=np.dtype(dtype_str), npairs=npairs)


def _ham_chunk(rng):
    i0, i1 = rng
    M, n, dtype = _HM["M"], _HM["n"], _HM["dtype"]
    shm = shared_memory.SharedMemory(name=_HM["shm_name"])
    out = np.ndarray((_HM["npairs"],), dtype=dtype, buffer=shm.buf)
    for i in range(i0, min(i1, n - 1)):
        B = M[i + 1:]                       # (n-i-1, R)
        diff = np.count_nonzero(B != M[i], axis=1).astype(dtype)
        off = _row_offset(i, n)
        out[off:off + diff.shape[0]] = diff
    shm.close()
    return None


def hamming_counts(M, jobs):
    """Pairwise disagreement counts (condensed, uint8/uint16) for rows of M."""
    n, R = M.shape
    npairs = n * (n - 1) // 2
    dtype = np.uint8 if R < 256 else np.uint16
    shm = shared_memory.SharedMemory(create=True, size=npairs * np.dtype(dtype).itemsize)
    try:
        # contiguous row-chunks; many chunks -> dynamic load balancing
        chunk = max(64, n // (jobs * 8))
        ranges = [(i, min(i + chunk, n)) for i in range(0, n, chunk)]
        if jobs <= 1:
            _init_ham(M, shm.name, n, R, dtype.__name__, npairs)
            for rng in ranges:
                _ham_chunk(rng)
        else:
            with ProcessPoolExecutor(max_workers=jobs, initializer=_init_ham,
                                     initargs=(M, shm.name, n, R, dtype.__name__, npairs)) as ex:
                list(ex.map(_ham_chunk, ranges))
        buf = np.ndarray((npairs,), dtype=dtype, buffer=shm.buf)
        return buf.copy(), dtype
    finally:
        shm.close()
        shm.unlink()


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run_cmmd(resolution_start=0.0, resolution_end=30.0, interval=0.5,
             seed=42, layer_weighting="normalized", nodelist=None,
             compute_distance=True, jobs=None, n_iterations=-1) -> Path:
    config.ensure_dirs()
    jobs = jobs or max(1, (os.cpu_count() or 2) // 2)

    layers, graphs, names = load_multiplex(config.HARMONIZED_DIR)
    n = len(names)
    print(f"[cmmd] {n:,} nodes across {len(layers)} layers: "
          + ", ".join(f"{l}({g.ecount():,}e)" for l, g in zip(layers, graphs)))
    print(f"[cmmd] weighting={layer_weighting}  jobs={jobs}")

    steps = int(round((resolution_end - resolution_start) / interval))
    resolutions = [round(resolution_start + i * interval, 10) for i in range(steps + 1)]
    print(f"[cmmd] {len(resolutions)} resolutions {resolutions[0]}..{resolutions[-1]} step {interval}")

    membership, communities = sweep(resolutions, config.HARMONIZED_DIR,
                                    layer_weighting, seed, jobs, n_iterations)
    print(f"[cmmd] communities per resolution: {communities[0]} .. {communities[-1]} "
          f"(min {min(communities)}, max {max(communities)})")

    # ---- membership outputs ----
    res_cols = [f"res_{g:g}" for g in resolutions]
    np.save(config.CMMD_DIR / "membership.npy", membership)
    (config.CMMD_DIR / "membership_nodes.txt").write_text("\n".join(names), encoding="utf-8")
    mem_path = config.CMMD_DIR / "membership_matrix.tsv"
    with open(mem_path, "w", newline="", encoding="utf-8") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["gene"] + res_cols + ["Pattern"])
        for i, name in enumerate(names):
            row = membership[i].tolist()
            w.writerow([name] + row + ["_".join(map(str, row))])
    print(f"[ok] membership_matrix.tsv + membership.npy ({n:,} x {len(resolutions)})")

    # ---- node selection ----
    if nodelist:
        wanted = {ln.strip() for ln in Path(nodelist).read_text(encoding="utf-8").splitlines() if ln.strip()}
        sel = [i for i, nm in enumerate(names) if nm in wanted]
        print(f"[cmmd] nodelist: {len(sel):,}/{len(wanted):,} genes found")
    else:
        sel = list(range(n))

    dist_info = {}
    if compute_distance:
        sel_names = [names[i] for i in sel]
        M = np.ascontiguousarray(membership[sel])
        pairs = len(sel) * (len(sel) - 1) // 2
        print(f"[cmmd] Hamming over {len(sel):,} genes ({pairs:,} pairs) on {jobs} workers ...")
        counts, dtype = hamming_counts(M, jobs)
        np.save(config.CMMD_DIR / "hamming_counts.npy", counts)
        (config.CMMD_DIR / "hamming_nodes.txt").write_text("\n".join(sel_names), encoding="utf-8")
        print(f"[ok] hamming_counts.npy ({counts.nbytes/1e6:.0f} MB, {dtype.__name__}) + hamming_nodes.txt")
        d = counts / len(resolutions)
        dist_info = {"genes": len(sel), "pairs": int(pairs), "dtype": dtype.__name__,
                     "n_resolutions": len(resolutions),
                     "distance_min": float(d.min()), "distance_max": float(d.max()),
                     "distance_mean": float(d.mean())}

    write_provenance(config.CMMD_DIR / "cmmd_report.json", {
        "params": {"resolution_start": resolution_start, "resolution_end": resolution_end,
                   "interval": interval, "seed": seed, "layer_weighting": layer_weighting,
                   "n_iterations": n_iterations, "nodelist": nodelist, "jobs": jobs},
        "layers": {l: g.ecount() for l, g in zip(layers, graphs)},
        "nodes": n,
        "n_resolutions": len(resolutions),
        "resolutions": resolutions,
        "communities_per_resolution": dict(zip([f"{g:g}" for g in resolutions], communities)),
        "distance": dist_info,
        "decode": "hamming_distance = hamming_counts / n_resolutions",
    })
    print(f"[ok] CmmD done -> {config.CMMD_DIR}")
    return config.CMMD_DIR


def main() -> None:
    ap = argparse.ArgumentParser(description="Run CmmD multilayer community detection.")
    ap.add_argument("--resolution-start", type=float, default=0.0)
    ap.add_argument("--resolution-end", type=float, default=30.0)
    ap.add_argument("--interval", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--layer-weighting", choices=["normalized", "equal"], default="normalized")
    ap.add_argument("--n-iterations", type=int, default=-1,
                    help="Leiden iterations per resolution (-1 = to convergence; 2 = fast)")
    ap.add_argument("--nodelist", default=None, help="file of gene symbols to restrict the distance matrix")
    ap.add_argument("--jobs", type=int, default=None, help="parallel worker processes (default: half of logical CPUs)")
    ap.add_argument("--no-distance", action="store_true", help="skip the Hamming distance matrix")
    args = ap.parse_args()
    run_cmmd(resolution_start=args.resolution_start, resolution_end=args.resolution_end,
             interval=args.interval, seed=args.seed, layer_weighting=args.layer_weighting,
             nodelist=args.nodelist, compute_distance=not args.no_distance, jobs=args.jobs,
             n_iterations=args.n_iterations)


if __name__ == "__main__":
    main()
