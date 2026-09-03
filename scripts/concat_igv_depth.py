#!/usr/bin/env python3
"""
Assemble one combined image per variant: an SV-info header banner on top, then the IGV
screenshot, then (for DEL/DUP) the long-read depth plot. Plots are paired by their shared
basename `{family}_{ID}.png`, which both the IGV track (makeigvpesr.py) and the depth track
(plot_longread_depth.py) emit.

Variants without a depth plot (SNV/indel, non-DEL/DUP SVs) carry their IGV plot through with
just the header, so every plotted variant gets the same header (ID, length, SVTYPE,
coords). The header text is looked up from --varfile (the canonical 6-col BED) by matching
each image's ID; if no varfile is given, images are emitted without a header.
"""

import argparse
import glob
import gzip
import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

import draw_pedigree
import sample_labels

HEADER_H = 60
PED_H = 500           # target pedigree height in the header band
HEADER_PAD = 16
_FONT_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
]


def _load_font(size=28):
    for path in _FONT_CANDIDATES:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def trim_bottom_bg(img, pad=6, white=240, line_frac=0.6):
    """Crop the empty region IGV leaves below the last track. IGV panels sit on a near-white
    background (~250) with a slightly different-shaded left name gutter; the empty area below
    the last track holds only thin vertical dividers/markers (the gutter border, region and
    center lines) that run its full height. So: treat any near-white pixel as background,
    drop columns that are non-background for >= `line_frac` of the HEIGHT (those vertical
    lines), then crop just below the lowest remaining real-content row (reads, coverage,
    gene models, a track label, or a track's horizontal separator).

    Height/track-count agnostic: more family members add read tracks higher up but the empty
    tail is detected the same way, and the reliable vertical line (the gutter border, which
    spans H minus the ~130px ruler) only becomes a LARGER height fraction as H grows.
    Dropping a thin column never empties a horizontally-spanning content row, so this can
    under-trim (harmless) but won't clip real tracks. Width and top are untouched."""
    arr = np.asarray(img.convert("RGB"))
    H, W, _ = arr.shape
    nonbg = arr.min(axis=2) < white
    col_is_line = nonbg.sum(axis=0) >= line_frac * H     # thin full-height dividers/markers
    content = nonbg.copy()
    content[:, col_is_line] = False
    nz = np.nonzero(np.any(content, axis=1))[0]
    if nz.size == 0:
        return img
    last = min(H, int(nz[-1]) + 1 + pad)
    return img if last >= H else img.crop((0, 0, W, last))


def vstack_imgs(top, bottom):
    # align to the WIDEST pane (upscale the narrower one) so the high-res IGV snapshot is
    # never shrunk down to the depth plot's width
    w = max(top.width, bottom.width)
    scale = lambda im: im if im.width == w else im.resize((w, int(im.height * w / im.width)))
    top, bottom = scale(top), scale(bottom)
    out = Image.new("RGB", (w, top.height + bottom.height), (255, 255, 255))
    out.paste(top, (0, 0))
    out.paste(bottom, (0, top.height))
    return out


def hstack_imgs(left, right):
    # align to the TALLEST pane (upscale the shorter one) so neither snapshot is shrunk
    h = max(left.height, right.height)
    scale = lambda im: im if im.height == h else im.resize((int(im.width * h / im.height), h))
    left, right = scale(left), scale(right)
    out = Image.new("RGB", (left.width + right.width, h), (255, 255, 255))
    out.paste(left, (0, 0))
    out.paste(right, (left.width, 0))
    return out


def load_variant_info(varfile):
    """{ID: (chrom, start, end, svtype, carriers)} from the canonical 6-col BED (maybe gzipped);
    carriers is the set of sample IDs in col6 (used to locate the variant's family)."""
    info = {}
    opener = gzip.open if varfile.endswith(".gz") else open
    with opener(varfile, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                carriers = set(f[5].split(",")) if len(f) > 5 and f[5] else set()
                info[f[3]] = (f[0], int(f[1]), int(f[2]), f[4], carriers)
            except ValueError:
                continue
    return info


def load_genotypes(path):
    """{ID: {sample: GT}} from the reconciled 'ID<TAB>sample=GT,...' table (maybe gzipped)."""
    gts = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            d = {}
            if len(f) > 1 and f[1]:
                for tok in f[1].split(","):
                    tok = tok.strip()
                    if "=" in tok:
                        s, gt = tok.split("=", 1)
                        d[s.strip()] = gt.strip()
            gts[f[0]] = d
    return gts


def build_ped_context(ped_path, genotypes_path):
    """Bundle everything the pedigree glyph needs: family rows, per-family role labels, a
    sample->family index, and the per-variant genotype map. Returns None if no ped is given."""
    if not ped_path:
        return None
    fams = sample_labels.read_ped_families(ped_path)
    roles = sample_labels.assign_labels(ped_path)
    sample_family = {m["iid"]: fam for fam, members in fams.items() for m in members}
    genotypes = load_genotypes(genotypes_path) if genotypes_path else {}
    return {"fams": fams, "roles": roles,
            "sample_family": sample_family, "genotypes": genotypes}


def variant_id_for(stem, info):
    """The variant ID embedded in a `{family}_{ID}` stem: the stem itself if it is an ID,
    else the longest known ID that is a trailing `_ID` of it. None if unmatched."""
    if stem in info:
        return stem
    cands = [i for i in info if stem.endswith("_" + i)]
    return max(cands, key=len) if cands else None


def reordered_name(stem, info):
    """`{family}_{ID}` -> `{ID}_{family}.png` so combined outputs sort by variant ID.
    Falls back to the original stem if the ID can't be identified."""
    vid = variant_id_for(stem, info)
    if not vid or stem == vid:
        return stem + ".png"
    fam = stem[:-(len(vid) + 1)]        # strip trailing "_<vid>"
    return f"{vid}_{fam}.png" if fam else stem + ".png"


def header_text_for(stem, info):
    """Format the SV-info header for a stem, or None if its variant ID is unmatched."""
    vid = variant_id_for(stem, info)
    if vid is None:
        return None
    chrom, start, end, svtype = info[vid][:4]
    return f"{vid}   {end - start:,} bp   {svtype}   {chrom}:{start}-{end}"


def pedigree_for(stem, info, ped_ctx):
    """Pedigree image for a variant's family (located via its carriers), annotated with the
    variant's genotypes + carriers, or None if there is no ped or no resolvable family."""
    if not ped_ctx:
        return None
    vid = variant_id_for(stem, info)
    if vid is None:
        return None
    carriers = info[vid][4]
    fam = next((ped_ctx["sample_family"][c] for c in sorted(carriers)
                if c in ped_ctx["sample_family"]), None)
    if fam is None:
        return None
    return draw_pedigree.draw_pedigree(
        ped_ctx["fams"].get(fam, []), ped_ctx["roles"].get(fam, {}),
        gts=ped_ctx["genotypes"].get(vid, {}), carriers=carriers, target_h=PED_H)


def build_header(stem, info, width, ped_ctx):
    """A header band (width `width`): SV-info text on the left, pedigree on the right. Returns
    None when there is neither text nor pedigree to show."""
    text = header_text_for(stem, info) if info else None
    ped_img = pedigree_for(stem, info, ped_ctx)
    if text is None and ped_img is None:
        return None
    if ped_img is not None and ped_img.width > width * 0.48:      # keep room for the text
        max_w = max(1, int(width * 0.48))
        ped_img = ped_img.resize((max_w, max(1, int(ped_img.height * max_w / ped_img.width))))

    height = max(HEADER_H, (ped_img.height + HEADER_PAD) if ped_img is not None else 0)
    header = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(header)
    if text:
        draw.text((20, height // 2 - 20), text, fill=(0, 0, 0), font=_load_font(40))
    if ped_img is not None:
        header.paste(ped_img, (max(0, width - ped_img.width - HEADER_PAD),
                               max(0, (height - ped_img.height) // 2)))
    return header


def finalize(body, stem, info, out_path, ped_ctx=None):
    header = build_header(stem, info, body.width, ped_ctx)
    if header is not None:
        body = vstack_imgs(header, body)
    body.save(out_path)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--igv-dir", required=True, help="dir of IGV {family}_{ID}.png plots")
    ap.add_argument("--depth-dir", required=True, help="dir of depth {family}_{ID}.png plots")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--varfile", default=None,
                    help="canonical 6-col BED for the SV-info header (chrom,start,end,ID,svtype,samples)")
    ap.add_argument("--ped", default=None,
                    help="ped file; enables the per-variant pedigree glyph in the header")
    ap.add_argument("--genotypes", default=None,
                    help="reconciled 'ID<TAB>sample=GT,...' table to annotate the pedigree")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    info = load_variant_info(args.varfile) if args.varfile else {}
    ped_ctx = build_ped_context(args.ped, args.genotypes)

    igv_plots = collect_igv_plots(args.igv_dir)
    depth_plots = {os.path.basename(p): p for p in glob.glob(os.path.join(args.depth_dir, "*.png"))}

    n_combined = n_igv_only = 0
    for name, igv_path in sorted(igv_plots.items()):
        stem = name[:-len(".png")]
        out_path = os.path.join(args.outdir, reordered_name(stem, info))
        igv_img = trim_bottom_bg(Image.open(igv_path).convert("RGB"))
        if name in depth_plots:
            body = vstack_imgs(igv_img, Image.open(depth_plots[name]).convert("RGB"))
            n_combined += 1
        else:
            body = igv_img
            n_igv_only += 1
        finalize(body, stem, info, out_path, ped_ctx)

    # depth plots with no matching IGV plot (shouldn't happen, but don't drop them silently)
    for name, depth_path in sorted(depth_plots.items()):
        if name not in igv_plots:
            stem = name[:-len(".png")]
            finalize(Image.open(depth_path).convert("RGB"), stem, info,
                     os.path.join(args.outdir, reordered_name(stem, info)), ped_ctx)
            sys.stderr.write(f"WARNING: depth plot {name} had no matching IGV plot\n")

    sys.stderr.write(f"Combined {n_combined}, IGV-only {n_igv_only}, written to {args.outdir}\n")


def collect_igv_plots(igv_dir):
    """Return {`{stem}.png`: path}, first merging any `{stem}.left.png` + `{stem}.right.png`
    breakpoint panes (emitted by makeigvpesr.py for large SVs) into one `{stem}.png` so they
    can pair with the single depth plot of the same name. Left/right are placed SIDE BY SIDE
    (left | right) at full resolution; the depth track below is upscaled to span both."""
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
        imgs = [Image.open(p).convert("RGB") for p in ordered]
        merged = imgs[0] if len(imgs) == 1 else hstack_imgs(imgs[0], imgs[1])
        merged.save(out)
        plain[stem] = out

    return {stem + ".png": path for stem, path in plain.items()}


if __name__ == "__main__":
    main()
