#!/usr/bin/env python3
"""
Long-read read-depth plot for CNV (DEL/DUP) loci — the replacement for GATK-SV's RdTest
in the long-read visualization workflow.

For each variant in the per-family BED it draws a normalized-depth profile across the
region +/- flank, one line per family sample (carriers highlighted, others gray), shades
the called interval, and marks the expected het DEL (0.5x) / DUP (1.5x) levels. One PNG per
variant is written to --outdir, named {family}_{ID}.png so it pairs with the IGV plot of
the same variant.

Depth comes from mosdepth `--by <windows>` output: one {sample}.regions.bed.gz per sample,
columns chrom,start,end,depth (mean depth per window). Normalization is local and per
sample: divide by the median depth over that sample's windows lying in the flanks (outside
the called interval), which needs no genome-wide coverage stats.
"""

import argparse
import glob
import gzip
import os
import sys
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_ped(ped_file):
    """Return dict sample -> {'family','affected'} using the first 6 standard columns."""
    info = {}
    with open(ped_file) as fh:
        for line in fh:
            if line.startswith(("#", "FamilyID", "family_id")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            fam, iid, _mo, _fa, _sex, aff = f[:6]
            info[iid] = {"family": fam, "affected": aff}
    return info


def load_depth(path):
    """Load a mosdepth regions.bed.gz -> list of (chrom, start, end, depth)."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, d = line.rstrip("\n").split("\t")[:4]
            rows.append((c, int(s), int(e), float(d)))
    return rows


def sample_from_filename(path):
    base = os.path.basename(path)
    return base.split(".regions.bed.gz")[0]


def norm_factor(windows, chrom, cnv_start, cnv_end):
    """Median depth over flank windows (outside the called interval) on this chrom."""
    flank = [d for (c, s, e, d) in windows
             if c == chrom and (e <= cnv_start or s >= cnv_end)]
    pool = flank if len(flank) >= 3 else [d for (c, s, e, d) in windows if c == chrom]
    pool = [d for d in pool if d > 0]
    return median(pool) if pool else 1.0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True,
                    help="per-family variant BED: chrom,start,end,ID,svtype,samples")
    ap.add_argument("--ped", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--flank", type=int, default=5000,
                    help="bp of flank shown on each side of the interval")
    ap.add_argument("--depth-dir", default=".",
                    help="directory containing {sample}.regions.bed.gz files")
    ap.add_argument("--outdir", default="rd_plots")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ped = read_ped(args.ped)

    depth = {}
    for path in glob.glob(os.path.join(args.depth_dir, "*.regions.bed.gz")):
        depth[sample_from_filename(path)] = load_depth(path)
    if not depth:
        sys.exit(f"ERROR: no *.regions.bed.gz found in {args.depth_dir}")

    n_plots = 0
    with open(args.bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end, vid, svtype = f[0], int(f[1]), int(f[2]), f[3], f[4].upper()
            carriers = set(f[5].split(",")) if len(f) > 5 and f[5] else set()
            if svtype not in ("DEL", "DUP"):
                continue

            lo, hi = max(0, start - args.flank), end + args.flank
            fig, ax = plt.subplots(figsize=(10, 4))

            for sample, windows in sorted(depth.items()):
                factor = norm_factor(windows, chrom, start, end)
                xs, ys = [], []
                for (c, s, e, d) in windows:
                    if c == chrom and e > lo and s < hi:
                        xs.append((s + e) / 2.0)
                        ys.append(d / factor if factor else 0.0)
                if not xs:
                    continue
                pts = sorted(zip(xs, ys))
                xs, ys = [p[0] for p in pts], [p[1] for p in pts]
                is_carrier = sample in carriers
                ax.plot(xs, ys,
                        color="crimson" if is_carrier else "0.7",
                        lw=2.0 if is_carrier else 1.0,
                        zorder=3 if is_carrier else 1,
                        label=f"{sample} (carrier)" if is_carrier else sample)

            ax.axvspan(start, end, color="lightblue", alpha=0.3, zorder=0)
            ax.axhline(1.0, color="black", ls="--", lw=0.8, zorder=2)
            ax.axhline(0.5 if svtype == "DEL" else 1.5, color="green", ls=":", lw=0.8, zorder=2)
            ax.set_ylim(0, 3)
            ax.set_xlabel(f"{chrom} position")
            ax.set_ylabel("normalized depth")
            ax.set_title(f"{vid}  {svtype}  {chrom}:{start}-{end}")
            ax.legend(fontsize=6, ncol=2, loc="upper right")
            fig.tight_layout()
            fig.savefig(os.path.join(args.outdir, f"{args.family}_{vid}.png"), dpi=120)
            plt.close(fig)
            n_plots += 1

    sys.stderr.write(f"Wrote {n_plots} depth plots to {args.outdir}\n")


if __name__ == "__main__":
    main()
