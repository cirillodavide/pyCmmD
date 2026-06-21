"""Lazy accessor for the CmmD pairwise Hamming distance matrix.

The distance is stored compressed as condensed uint8 disagreement counts; this
class memory-maps it (so nothing big is loaded into RAM) and decodes on access:
distance = count / n_resolutions.

    d = CmmDDistance()                 # loads from data/cmmd/
    d.dist("BRCA1", "BARD1")           # single pair, O(1)
    d.row("TP53")                      # distances to all genes
    d.nearest("TP53", k=10)           # k most persistently co-clustered genes

To materialize the full (decoded) distance matrix from the compressed counts:

    python -m pycmmd.distance --form square          # -> data/cmmd/hamming_square.npy
    python -m pycmmd.distance --form condensed       # -> hamming_distance_condensed.npy
    python -m pycmmd.distance --form square --tsv     # also a labeled TSV (large/slow)

Row/column labels are in hamming_nodes.txt (same order as the matrix).
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from . import config


class CmmDDistance:
    def __init__(self, cmmd_dir=None):
        d = cmmd_dir or config.CMMD_DIR
        self.nodes = (d / "hamming_nodes.txt").read_text(encoding="utf-8").splitlines()
        self.index = {g: i for i, g in enumerate(self.nodes)}
        self.n = len(self.nodes)
        self.counts = np.load(d / "hamming_counts.npy", mmap_mode="r")
        self.R = json.loads((d / "cmmd_report.json").read_text(encoding="utf-8"))["n_resolutions"]

    def _pos(self, i: int, j: int) -> int:
        if i > j:
            i, j = j, i
        n = self.n
        return n * (n - 1) // 2 - (n - i) * (n - i - 1) // 2 + (j - i - 1)

    def dist(self, a: str, b: str) -> float:
        i, j = self.index[a], self.index[b]
        if i == j:
            return 0.0
        return int(self.counts[self._pos(i, j)]) / self.R

    def row(self, gene: str) -> np.ndarray:
        """Distances from `gene` to every node (self = 0), as float32."""
        i = self.index[gene]
        n = self.n
        out = np.empty(n, dtype=np.float32)
        out[i] = 0.0
        # pairs (j, i) for j < i, and (i, j) for j > i
        for j in range(n):
            if j == i:
                continue
            out[j] = int(self.counts[self._pos(i, j)])
        out /= self.R
        return out

    def nearest(self, gene: str, k: int = 10):
        """Return [(gene, distance), ...] of the k closest genes (excluding self)."""
        d = self.row(gene)
        i = self.index[gene]
        order = np.argsort(d, kind="stable")
        out = [(self.nodes[j], float(d[j])) for j in order if j != i]
        return out[:k]

    # -- full-matrix materialization ---------------------------------------- #
    def to_condensed(self, dtype=np.float32) -> np.ndarray:
        """Decoded condensed distance vector (len n*(n-1)/2)."""
        return np.asarray(self.counts).astype(dtype) / self.R

    def to_square(self, dtype=np.float32) -> np.ndarray:
        """Decoded full n x n distance matrix (zero diagonal)."""
        from scipy.spatial.distance import squareform
        sq = squareform(np.asarray(self.counts))   # uint8 n x n, diag 0
        return sq.astype(dtype) / self.R

    def to_tsv(self, path, dtype=np.float32) -> None:
        """Stream a labeled n x n distance TSV (genes as header + row labels)."""
        n = self.n
        with open(path, "w", encoding="utf-8") as fo:
            fo.write("gene\t" + "\t".join(self.nodes) + "\n")
            for i in range(n):
                rowvals = self.row(self.nodes[i]).astype(dtype)
                fo.write(self.nodes[i] + "\t" + "\t".join(f"{v:g}" for v in rowvals) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Materialize the full CmmD distance matrix from compressed counts.")
    ap.add_argument("--form", choices=["square", "condensed"], default="square")
    ap.add_argument("--dtype", default="float32", help="float32 (default) or float64")
    ap.add_argument("--out", default=None, help="output .npy path (default in data/cmmd/)")
    ap.add_argument("--tsv", action="store_true",
                    help="also write a labeled TSV (square only; large and slow)")
    args = ap.parse_args()

    d = CmmDDistance()
    dtype = np.dtype(args.dtype)
    n, R = d.n, d.R
    if args.form == "square":
        gib = n * n * dtype.itemsize / 1e9
        print(f"[distance] building {n} x {n} square ({dtype.name}, ~{gib:.1f} GB in RAM) ...")
        arr = d.to_square(dtype)
        out = args.out or (config.CMMD_DIR / "hamming_square.npy")
    else:
        gib = (n * (n - 1) // 2) * dtype.itemsize / 1e9
        print(f"[distance] building condensed vector ({dtype.name}, ~{gib:.1f} GB) ...")
        arr = d.to_condensed(dtype)
        out = args.out or (config.CMMD_DIR / "hamming_distance_condensed.npy")
    np.save(out, arr)
    print(f"[ok] saved {out}  shape={arr.shape} dtype={arr.dtype} "
          f"({arr.nbytes/1e6:.0f} MB on disk)  decode=counts/{R}")

    if args.tsv and args.form == "square":
        tsv = str(out).replace(".npy", ".tsv")
        print(f"[distance] writing labeled TSV {tsv} (this is large/slow) ...")
        d.to_tsv(tsv, dtype)
        print(f"[ok] saved {tsv}")


if __name__ == "__main__":
    main()
