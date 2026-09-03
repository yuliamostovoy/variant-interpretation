#!/usr/bin/env python3
"""
Draw a compact pedigree for one family as a PIL image, for the combined-plot header.

Standard pedigree conventions: male = square, female = circle, unknown sex = diamond; a filled
symbol = affected, open = unaffected; a crimson outline marks a carrier of the variant. Each
individual is labeled below its symbol with its relationship role (from sample_labels.py) + ID
and, when a genotype map is supplied, its raw VCF genotype on the next line.

`draw_pedigree(family_rows, labels, gts=None, carriers=None, target_h=...)` returns a
`PIL.Image`. Layout is generational (founders on top, children below their parents, with mating
and sibship connectors) and scoped to the common founder-couple + children shapes; anything that
doesn't resolve cleanly falls back to a simple one-row-per-generation grid, and any failure at
all degrades to a blank image so the caller never crashes.

family_rows: list of dicts with keys iid, fa, mo, sex, aff (raw ped strings; missing parent = 0).
"""

import argparse
import io
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, RegularPolygon
from PIL import Image

_MISSING = {"0", "", ".", "-9"}

R = 0.22            # symbol half-size (data units)
DY = 1.5           # vertical spacing between generations
DX = 1.4           # horizontal spacing between individuals
DPI = 150          # native render resolution (kept high so the header glyph stays crisp)
FONTSIZE = 16      # node label / genotype font (points)


def _missing(x):
    return x is None or x in _MISSING


def _generations(members, by_iid):
    """iid -> generation index (founders = 0), by iterative propagation. Raises on cycles."""
    def parents(m):
        return [p for p in (m["fa"], m["mo"]) if not _missing(p) and p in by_iid]

    gen = {m["iid"]: 0 for m in members}
    for _ in range(len(members) + 1):
        changed = False
        for m in members:
            ps = parents(m)
            if ps:
                g = 1 + max(gen[p] for p in ps)
                if g != gen[m["iid"]]:
                    gen[m["iid"]] = g
                    changed = True
        if not changed:
            break
    else:
        raise ValueError("pedigree generation propagation did not converge (cycle?)")
    return gen


def _layout(members, by_iid, gen):
    """iid -> (x, y). Couples (co-parents of a shared child) are placed adjacent and centered
    over their children; children are centered under their parents. Bottom-up."""
    children = {}
    for m in members:
        for p in (m["fa"], m["mo"]):
            if not _missing(p) and p in by_iid:
                children.setdefault(p, []).append(m["iid"])

    # co-parent of each individual (other parent of any shared child), if present in the family
    spouse = {}
    for m in members:
        ps = [p for p in (m["fa"], m["mo"]) if not _missing(p) and p in by_iid]
        if len(ps) == 2:
            spouse.setdefault(ps[0], set()).add(ps[1])
            spouse.setdefault(ps[1], set()).add(ps[0])

    # a married-in spouse is drawn on their mate's generation, not their own founder row
    gen = dict(gen)
    for _ in range(len(members) + 1):
        changed = False
        for a, mates in spouse.items():
            for b in mates:
                g = max(gen[a], gen[b])
                if gen[a] != g:
                    gen[a] = g
                    changed = True
                if gen[b] != g:
                    gen[b] = g
                    changed = True
        if not changed:
            break

    order = {m["iid"]: i for i, m in enumerate(members)}
    by_gen = {}
    for iid, g in gen.items():
        by_gen.setdefault(g, []).append(iid)

    xpos = {}
    for g in sorted(by_gen, reverse=True):                 # deepest generation first
        nodes = sorted(by_gen[g], key=lambda i: order[i])
        used = set()
        units = []                                         # each unit is a list of iids
        for iid in nodes:
            if iid in used:
                continue
            mates = [s for s in spouse.get(iid, ()) if s in gen and gen[s] == g]
            if mates:
                mate = min(mates, key=lambda i: order[i])
                pair = sorted([iid, mate], key=lambda i: order[i])
                units.append(pair)
                used.update(pair)
            else:
                units.append([iid])
                used.add(iid)

        def desired(unit):
            kids = []
            for member in unit:
                kids += [k for k in children.get(member, []) if k in xpos]
            return sum(xpos[k] for k in kids) / len(kids) if kids else None

        anchored = [u for u in units if desired(u) is not None]
        free = [u for u in units if desired(u) is None]
        anchored.sort(key=desired)

        cursor = None
        for u in anchored:
            w = (len(u) - 1) * DX
            center = desired(u)
            if cursor is not None:
                center = max(center, cursor + DX + w / 2.0)
            start = center - w / 2.0
            for i, member in enumerate(u):
                xpos[member] = start + i * DX
            cursor = center + w / 2.0
        for u in free:
            start = 0.0 if cursor is None else cursor + DX
            for i, member in enumerate(u):
                xpos[member] = start + i * DX
            cursor = start + (len(u) - 1) * DX

    return {iid: (xpos[iid], -gen[iid] * DY) for iid in xpos}, children


def _grid_layout(members, gen):
    """Fallback: one row per generation, individuals spread left-to-right in ped order."""
    by_gen = {}
    for i, m in enumerate(members):
        by_gen.setdefault(gen.get(m["iid"], 0), []).append((i, m["iid"]))
    pos = {}
    for g, row in by_gen.items():
        for col, (_, iid) in enumerate(sorted(row)):
            pos[iid] = (col * DX, -g * DY)
    return pos


def _draw_symbol(ax, x, y, sex, affected, carrier):
    face = "black" if affected else "white"
    edge = "crimson" if carrier else "black"
    lw = 2.6 if carrier else 1.2
    if sex == "1":                                         # male: square
        ax.add_patch(Rectangle((x - R, y - R), 2 * R, 2 * R,
                               facecolor=face, edgecolor=edge, linewidth=lw, zorder=3))
    elif sex == "2":                                       # female: circle
        ax.add_patch(Circle((x, y), R, facecolor=face, edgecolor=edge,
                            linewidth=lw, zorder=3))
    else:                                                  # unknown: diamond
        ax.add_patch(RegularPolygon((x, y), numVertices=4, radius=R * 1.35,
                                    orientation=0, facecolor=face, edgecolor=edge,
                                    linewidth=lw, zorder=3))


def _fig_to_image(fig, target_h):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI, bbox_inches="tight", pad_inches=0.05,
                facecolor="white")
    plt.close(fig)
    buf.seek(0)
    img = Image.open(buf).convert("RGB")
    if target_h and img.height != target_h:
        w = max(1, int(round(img.width * target_h / img.height)))
        img = img.resize((w, target_h))
    return img


def _blank(target_h):
    return Image.new("RGB", (max(1, target_h), max(1, target_h or 1)), (255, 255, 255))


def draw_pedigree(family_rows, labels, gts=None, carriers=None, target_h=240):
    """Render one family's pedigree to a PIL.Image; never raises."""
    gts = gts or {}
    carriers = set(carriers or ())
    labels = labels or {}
    try:
        members = [dict(m) for m in family_rows if m.get("iid")]
        if not members:
            return _blank(target_h)
        by_iid = {m["iid"]: m for m in members}

        try:
            gen = _generations(members, by_iid)
            pos, children = _layout(members, by_iid, gen)
            draw_links = True
        except Exception:
            gen = {m["iid"]: 0 for m in members}
            pos = _grid_layout(members, gen)
            children, draw_links = {}, False

        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        ax.axis("off")

        if draw_links:
            _draw_connectors(ax, members, by_iid, pos, children)

        for m in members:
            iid = m["iid"]
            if iid not in pos:
                continue
            x, y = pos[iid]
            _draw_symbol(ax, x, y, m.get("sex", "0"), m.get("aff") == "2", iid in carriers)
            # Label each node with the SHORT role (Pro/Mo/Fa/...) only -- the full sample ID is
            # long and already shown on the IGV track and depth legend, which the role maps to;
            # printing it here just collides with the neighbouring nodes. Fall back to a
            # truncated ID only when a sample has no role.
            role = labels.get(iid, "")
            label = role or (iid[:10] + "..." if len(iid) > 12 else iid)
            gt = gts.get(iid)
            text = label + ("\n" + gt if gt else "")
            ax.text(x, y - R - 0.06, text, ha="center", va="top",
                    fontsize=FONTSIZE, fontweight="bold", linespacing=1.2, zorder=4)

        xs = [p[0] for p in pos.values()]
        ys = [p[1] for p in pos.values()]
        ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
        ax.set_ylim(min(ys) - DY, max(ys) + R + 0.4)
        return _fig_to_image(fig, target_h)
    except Exception as e:
        sys.stderr.write("WARNING: pedigree render failed: {}\n".format(e))
        return _blank(target_h)


def _draw_connectors(ax, members, by_iid, pos, children):
    line = dict(color="black", linewidth=1.0, zorder=1)
    drawn_couples = set()
    for m in members:
        iid = m["iid"]
        kids = [k for k in children.get(iid, []) if k in pos]
        if not kids:
            continue
        # group this individual's children by the co-parent so each couple draws once
        by_mate = {}
        for k in kids:
            km = by_iid[k]
            mates = [p for p in (km["fa"], km["mo"])
                     if not _missing(p) and p in by_iid and p != iid]
            by_mate.setdefault(mates[0] if mates else None, []).append(k)
        for mate, mkids in by_mate.items():
            key = frozenset([iid, mate]) if mate else frozenset([iid, tuple(sorted(mkids))])
            if key in drawn_couples:
                continue
            drawn_couples.add(key)
            px, py = pos[iid]
            if mate and mate in pos:
                mx, my = pos[mate]
                ax.plot([px, mx], [py, my], **line)               # mating line
                mid_x, mid_y = (px + mx) / 2.0, (py + my) / 2.0
            else:
                mid_x, mid_y = px, py                              # single parent
            sib_y = mid_y - DY * 0.5
            ax.plot([mid_x, mid_x], [mid_y, sib_y], **line)       # drop to sibship line
            kxs = [pos[k][0] for k in mkids]
            ax.plot([min(kxs), max(kxs)], [sib_y, sib_y], **line)  # sibship line
            for k in mkids:
                kx, ky = pos[k]
                ax.plot([kx, kx], [sib_y, ky + R], **line)         # drop to each child


def _load_genotypes(path):
    gts = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            for tok in "\t".join(f[1:]).replace("\t", ",").split(","):
                tok = tok.strip()
                if "=" in tok:
                    s, gt = tok.split("=", 1)
                    gts[s.strip()] = gt.strip()
    return gts


def main():
    import sample_labels
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ped", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--carriers", default="", help="comma-separated carrier sample IDs")
    ap.add_argument("--genotypes", default=None,
                    help="TSV: any line 'x  sample=GT[,sample=GT...]' (ID column ignored)")
    ap.add_argument("--height", type=int, default=240)
    args = ap.parse_args()

    fams = sample_labels.read_ped_families(args.ped)
    labels = sample_labels.assign_labels(args.ped).get(args.family, {})
    rows = fams.get(args.family, [])
    carriers = {s for s in args.carriers.split(",") if s}
    gts = _load_genotypes(args.genotypes) if args.genotypes else {}
    img = draw_pedigree(rows, labels, gts=gts, carriers=carriers, target_h=args.height)
    img.save(args.out)
    sys.stderr.write("Wrote {} ({}x{})\n".format(args.out, img.width, img.height))


if __name__ == "__main__":
    main()
