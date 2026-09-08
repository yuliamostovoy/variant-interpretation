#!/usr/bin/env python3
"""
Lay out the long-read visualization combined plots for review in flipbook
(https://github.com/broadinstitute/flipbook).

flipbook shows every image in a subdirectory together on one page, so to get one review
page per variant we place each combined PNG in its own subfolder:

    <out>/<ID>_<family>/<ID>_<family>.png

and write a metadata table (<out>/flipbook_metadata.tsv) whose `Path` column is each variant's
relative subfolder (how flipbook keys metadata to a page) plus the SV facts, so each page
displays ID / type / coords / size / family / per-sample genotypes next to the image.

Usage:
    prepare_flipbook.py --combined-dir <dir of {ID}_{family}.png> \
        --varfile <canonical 6-col BED(.gz)> [--genotypes <ID<TAB>sample=GT,... (.gz)>] \
        --out <review dir>

Then review with (flipbook: https://github.com/broadinstitute/flipbook; `pip install flipbook`):
    cd <review dir>
    python3 -m flipbook -m flipbook_metadata.tsv -j /path/to/scripts/sv_review_form.json

The companion form schema scripts/sv_review_form.json defines the curation questions
(Call: Real/Artifact/Uncertain; Genotypes match evidence; Inheritance; free-text notes).
Responses are appended to flipbook_form_responses.tsv, keyed by each image's Path.
"""

import argparse
import gzip
import os
import shutil


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_varinfo(varfile):
    """{ID: (chrom, start, end, svtype, carriers_csv)} from the canonical 6-col BED."""
    info = {}
    with _open(varfile) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                int(f[1]); int(f[2])
            except ValueError:
                continue
            carriers = f[5] if len(f) > 5 else ""
            info[f[3]] = (f[0], int(f[1]), int(f[2]), f[4], carriers)
    return info


def load_gts(path):
    """{ID: {sample: GT}} parsed from the reconciled 'ID<TAB>sample=GT,...' table."""
    out = {}
    if not path:
        return out
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            d = {}
            if len(f) > 1 and f[1]:
                for tok in f[1].split(","):
                    if "=" in tok:
                        s, gt = tok.strip().split("=", 1)
                        d[s.strip()] = gt.strip()
            out[f[0]] = d
    return out


def load_family_members(ped_path):
    """{family: [iid, ...] in ped order} so genotypes can be shown for family members only."""
    fams = {}
    if not ped_path:
        return fams
    with open(ped_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "FamilyID", "family_id")):
                continue
            f = line.split()
            if len(f) < 2:
                continue
            fam, iid = f[0], f[1]
            members = fams.setdefault(fam, [])
            if iid not in members:
                members.append(iid)
    return fams


def render_gts(gt_map, members, carriers_csv):
    """Show GTs for the family's members (ped order) if known, else just the carriers."""
    ids = members or [c for c in carriers_csv.split(",") if c]
    return ", ".join(f"{s}={gt_map[s]}" for s in ids if s in gt_map)


def stem_to_id(stem, info):
    """The PNG is named {ID}_{family}.png; recover the ID by longest matching known ID prefix."""
    cands = [vid for vid in info if stem == vid or stem.startswith(vid + "_")]
    return max(cands, key=len) if cands else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combined-dir", required=True)
    ap.add_argument("--varfile", required=True)
    ap.add_argument("--genotypes", default=None)
    ap.add_argument("--ped", default=None,
                    help="restrict displayed genotypes to each variant's family members")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    info = load_varinfo(args.varfile)
    gts = load_gts(args.genotypes)
    fam_members = load_family_members(args.ped)

    os.makedirs(args.out, exist_ok=True)
    cols = ["Path", "ID", "SVTYPE", "size_bp", "locus", "family", "carriers", "genotypes"]
    rows = []
    pngs = sorted(f for f in os.listdir(args.combined_dir) if f.endswith(".png"))
    for png in pngs:
        stem = png[:-4]
        vid = stem_to_id(stem, info)
        family = stem[len(vid) + 1:] if vid and stem.startswith(vid + "_") else ""
        subdir = os.path.join(args.out, stem)
        os.makedirs(subdir, exist_ok=True)
        shutil.copy2(os.path.join(args.combined_dir, png), os.path.join(subdir, png))
        # flipbook keys the -m metadata table by each image's relative DIRECTORY (one row per
        # page), not the image file path -- so Path is the per-variant subfolder.
        rel = stem
        if vid and vid in info:
            chrom, start, end, svtype, carriers = info[vid]
            gt_str = render_gts(gts.get(vid, {}), fam_members.get(family, []), carriers)
            rows.append([rel, vid, svtype, str(end - start),
                         f"{chrom}:{start}-{end}", family, carriers, gt_str])
        else:
            rows.append([rel, stem, "", "", "", family, "", ""])

    with open(os.path.join(args.out, "flipbook_metadata.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    print(f"Prepared {len(rows)} variant page(s) under {args.out}")
    print(f"Metadata: {os.path.join(args.out, 'flipbook_metadata.tsv')}")


if __name__ == "__main__":
    main()
