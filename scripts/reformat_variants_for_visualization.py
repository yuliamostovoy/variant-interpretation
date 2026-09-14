#!/usr/bin/env python3
"""
Turn a candidate-variant list (generated from a VCF, typically via `bcftools query`) into
the canonical bgzipped 6-column BED the long-read visualization workflow consumes.

REQUIRED INPUT: a tab-separated file with these 6 columns, in this order:

    1. chrom    e.g. chr1
    2. pos      1-based start (use %POS in bcftools query)
    3. end      end coordinate (use %END in bcftools query)
    4. ID       variant ID (use %ID; "." is allowed and auto-filled)
    5. svtype   SV type (use %INFO/SVTYPE; "." allowed, becomes VAR — fine for SNV/indel)
    6. samples  comma-separated IDs of the samples to plot this variant for

A header line is optional (auto-detected and skipped). Example generator:

    bcftools view -i 'GT="alt"' calls.vcf.gz \\
      | bcftools query -f '%CHROM\\t%POS\\t%END\\t%ID\\t%INFO/SVTYPE\\t[%SAMPLE,]\\n' \\
      > candidates.tsv

(the trailing comma left by `[%SAMPLE,]` is stripped automatically.)

The `samples` IDs must match the individual IDs in the ped file and the sample->BAM map.
Output is plain TSV (chrom,start,end,ID,svtype,samples,svlen + header); the WDL task bgzips it.
The 7th column `svlen` is the variant's length in bp, taken from REF/ALT via --variant-info when
available (so insertions/deletions size by their allele length) and otherwise the coordinate
span. Consumers that want only the original 6 columns can ignore it.

OPTIONAL: pass --genotypes-raw to also emit a per-variant genotype table for the pedigree glyph.
It is a bcftools dump over the source VCF for ALL samples:

    bcftools query -f '%CHROM\\t%POS\\t%END\\t%ID\\t%INFO/SVTYPE[\\t%SAMPLE=%GT]\\n' calls.vcf.gz

Each row is reconciled to the SAME resolved variant ID written above -- by VCF ID when present,
falling back to locus (chrom,pos,end,svtype) for IDs synthesized from "." -- so the plots can join
genotypes by ID. The result (--genotypes-out, default <output>.genotypes.tsv) is
'ID<TAB>sample=GT,sample=GT,...'.
"""

import argparse
import sys


def is_header(fields):
    # header if the "start" column isn't an integer
    if len(fields) < 2:
        return True
    try:
        int(fields[1])
        return False
    except ValueError:
        return True


def clean_samples(raw):
    # accept comma/semicolon/space separators; drop empties; strip any "=GT" suffix
    for sep in (";", " "):
        raw = raw.replace(sep, ",")
    out = []
    for tok in raw.split(","):
        tok = tok.strip().split("=")[0]
        if tok:
            out.append(tok)
    return ",".join(out)


def norm_svtype(svtype):
    svtype = svtype.strip().upper()
    return "VAR" if svtype in ("", ".") else svtype


def resolve_id(chrom, start, end, vid, svtype):
    return vid if vid.strip() not in ("", ".") else f"{chrom}_{start}_{end}_{svtype}"


def allele_length(ref, alt):
    """Length in bp of a sequence-allele variant from its REF/ALT, or None for symbolic/
    breakend ALTs (the caller then uses the coordinate span). SNV/MNV -> len(REF); indel ->
    |len(REF) - len(ALT)| (the inserted or deleted bases, excluding the shared anchor base,
    e.g. G->GCA and GCA->G are both 2 bp)."""
    if not ref or not alt:
        return None
    alt = alt.split(",")[0].strip()      # first ALT if a record slipped through un-split
    ref = ref.strip()
    if not ref or not alt or alt == "." or alt[0] in "<*" or "[" in alt or "]" in alt:
        return None
    return len(ref) if len(ref) == len(alt) else abs(len(ref) - len(alt))


def load_allele_lengths(path):
    """From a VCF dump 'chrom, pos, end, ID, REF, ALT' (multiallelics pre-split) build
    (by_id, by_locus) lookups of allele length. The ID field may be ';'-joined for sites that
    were multiallelic, so register each sub-ID; the locus key is (chrom, pos, end)."""
    by_id, by_locus = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            length = allele_length(f[4], f[5])
            if length is None:
                continue
            try:
                by_locus.setdefault((f[0], int(f[1]), int(f[2])), length)
            except ValueError:
                continue
            for one in f[3].split(";"):
                one = one.strip()
                if one and one != ".":
                    by_id.setdefault(one, length)
    return by_id, by_locus


def write_genotypes(raw_path, out_path, id_by_locus, valid_ids):
    """Rewrite a bcftools GT dump (chrom, pos, end, ID, svtype, then one 'sample=GT' field per
    sample) to 'ID<TAB>sample=GT,...' keyed by the SAME resolved variant ID as the varfile, so
    the plots can join genotypes by ID. Match on the VCF ID first (robust to any residual
    coordinate-convention differences), then fall back to locus for variants whose ID was
    synthesized from ".". Variants not in the varfile are skipped."""
    n = 0
    with open(raw_path) as fh, open(out_path, "w") as out:
        out.write("#ID\tgenotypes\n")
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            rid = f[3].strip()
            if rid in valid_ids:
                vid = rid
            else:
                vid = id_by_locus.get((f[0], int(f[1]), int(f[2]), norm_svtype(f[4])))
            if vid is None:
                continue
            gts = ",".join(tok.strip() for tok in f[5:] if tok.strip())
            out.write(f"{vid}\t{gts}\n")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True,
                    help="6-column TSV (chrom,pos,end,ID,svtype,samples); header optional")
    ap.add_argument("--output", default="variants_for_visualization.bed",
                    help="output BED (plain TSV; bgzipped by the WDL task)")
    ap.add_argument("--variant-info", default=None,
                    help="optional bcftools dump 'chrom,pos,end,ID,REF,ALT' (multiallelics "
                         "split) used to size each variant by its allele length in the header")
    ap.add_argument("--genotypes-raw", default=None,
                    help="optional bcftools GT dump (chrom,pos,end,ID,svtype, then a "
                         "'sample=GT' field per sample) to reconcile against the resolved IDs")
    ap.add_argument("--genotypes-out", default=None,
                    help="where to write the reconciled 'ID<TAB>sample=GT,...' table "
                         "(default: <output>.genotypes.tsv)")
    args = ap.parse_args()

    rows = []
    with open(args.input) as fh:
        for lineno, line in enumerate(fh):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if lineno == 0 and is_header(f):
                continue
            if len(f) < 6:
                sys.exit(f"ERROR: line {lineno + 1} has {len(f)} columns, expected 6 "
                         f"(chrom,start,end,ID,svtype,samples): {line.rstrip()}")
            chrom, start, end, vid, svtype, samples = f[:6]
            start, end = int(start), int(end)
            samples = clean_samples(samples)
            svtype = norm_svtype(svtype)
            vid = resolve_id(chrom, start, end, vid, svtype)
            rows.append((chrom, start, end, vid, svtype, samples))

    rows.sort(key=lambda r: (r[0], r[1]))

    by_id, by_locus = load_allele_lengths(args.variant_info) if args.variant_info else ({}, {})

    with open(args.output, "w") as out:
        out.write("#chrom\tstart\tend\tID\tsvtype\tsamples\tsvlen\n")
        for chrom, start, end, vid, svtype, samples in rows:
            # allele length from REF/ALT (by ID, then locus); fall back to the coordinate span,
            # which is SVLEN for symbolic SVs. No flooring -- the value is the true length.
            svlen = by_id.get(vid)
            if svlen is None:
                svlen = by_locus.get((chrom, start, end))
            if svlen is None:
                svlen = end - start
            # An insertion occupies no reference span; callers still store END=POS+SVLEN, which
            # would off-center the plot and blow up the window. Collapse it to a point at POS
            # (svlen, computed above, keeps the true inserted length for the header).
            out_end = start if svtype == "INS" else end
            out.write(f"{chrom}\t{start}\t{out_end}\t{vid}\t{svtype}\t{samples}\t{svlen}\n")

    sys.stderr.write(f"Wrote {len(rows)} variants to {args.output}\n")

    if args.genotypes_raw:
        gt_out = args.genotypes_out or (args.output + ".genotypes.tsv")
        id_by_locus = {(c, s, e, sv): vid for c, s, e, vid, sv, _ in rows}
        valid_ids = {vid for _, _, _, vid, _, _ in rows}
        n = write_genotypes(args.genotypes_raw, gt_out, id_by_locus, valid_ids)
        sys.stderr.write(f"Wrote genotypes for {n} variants to {gt_out}\n")


if __name__ == "__main__":
    main()
