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


def hstack(paths, out_path):
    imgs = [Image.open(p).convert("RGB") for p in paths]
    h = min(im.height for im in imgs)
    scaled = [im if im.height == h else im.resize((int(im.width * h / im.height), h))
              for im in imgs]
    combined = Image.new("RGB", (sum(im.width for im in scaled), h), (255, 255, 255))
    x = 0
    for im in scaled:
        combined.paste(im, (x, 0))
        x += im.width
    combined.save(out_path)


def collect_igv_plots(igv_dir):
    """Return {`{stem}.png`: path}, first merging any `{stem}.left.png` + `{stem}.right.png`
    breakpoint panes (emitted by makeigvpesr.py for large SVs) into one `{stem}.png` so they
    can pair with the single depth plot of the same name."""
    plain = {}   # stem -> path
    sides = {}   # stem -> {"left": path, "right": path}
    for p in glob.glob(os.path.join(igv_dir, "*.png")):
        b = os.path.basename(p)
        if b.endswith(".left.png"):
            sides.setdefault(b[:-len(".left.png")], {})["left"] = p
        elif b.endswith(".right.png"):
            sides.setdefault(b[:-len(".right.png")], {})["right"] = p
        else:
            plain[b[:-len(".png")]] = p

    for stem, d in sides.items():
        if stem in plain:            # a single-pane plot already exists; prefer it
            continue
        ordered = [d[k] for k in ("left", "right") if k in d]
        out = os.path.join(igv_dir, stem + ".png")
        if len(ordered) == 1:
            Image.open(ordered[0]).convert("RGB").save(out)
        else:
            hstack(ordered, out)     # left | right
        plain[stem] = out

    return {stem + ".png": path for stem, path in plain.items()}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--igv-dir", required=True, help="dir of IGV {family}_{ID}.png plots")
    ap.add_argument("--depth-dir", required=True, help="dir of depth {family}_{ID}.png plots")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    igv_plots = collect_igv_plots(args.igv_dir)
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
