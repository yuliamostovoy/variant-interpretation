#!/usr/bin/env python3
"""
Stack the IGV screenshot (top) and the long-read depth plot (bottom) for each variant that
has both, into a single combined image. Pairs plots by their shared basename
`{family}_{ID}.png`, which both the IGV track (makeigvpesr.py) and the depth track
(plot_longread_depth.py) emit.

Variants without a depth plot (e.g. SNV/indel, SV that isn't DEL/DUP) simply carry their
IGV plot through unchanged, so the output directory holds one image per plotted variant.
"""

import argparse
import glob
import os
import sys

from PIL import Image


def vstack(top_path, bottom_path, out_path):
    top = Image.open(top_path).convert("RGB")
    bottom = Image.open(bottom_path).convert("RGB")
    w = min(top.width, bottom.width)

    def scale(img):
        if img.width == w:
            return img
        h = int(img.height * w / img.width)
        return img.resize((w, h))

    top, bottom = scale(top), scale(bottom)
    combined = Image.new("RGB", (w, top.height + bottom.height), (255, 255, 255))
    combined.paste(top, (0, 0))
    combined.paste(bottom, (0, top.height))
    combined.save(out_path)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--igv-dir", required=True, help="dir of IGV {family}_{ID}.png plots")
    ap.add_argument("--depth-dir", required=True, help="dir of depth {family}_{ID}.png plots")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    igv_plots = {os.path.basename(p): p for p in glob.glob(os.path.join(args.igv_dir, "*.png"))}
    depth_plots = {os.path.basename(p): p for p in glob.glob(os.path.join(args.depth_dir, "*.png"))}

    n_combined = n_igv_only = 0
    for name, igv_path in sorted(igv_plots.items()):
        out_path = os.path.join(args.outdir, name)
        if name in depth_plots:
            vstack(igv_path, depth_plots[name], out_path)
            n_combined += 1
        else:
            Image.open(igv_path).convert("RGB").save(out_path)
            n_igv_only += 1

    # depth plots with no matching IGV plot (shouldn't happen, but don't drop them silently)
    for name, depth_path in sorted(depth_plots.items()):
        if name not in igv_plots:
            Image.open(depth_path).convert("RGB").save(os.path.join(args.outdir, name))
            sys.stderr.write(f"WARNING: depth plot {name} had no matching IGV plot\n")

    sys.stderr.write(f"Combined {n_combined}, IGV-only {n_igv_only}, written to {args.outdir}\n")


if __name__ == "__main__":
    main()
