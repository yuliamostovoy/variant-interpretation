#!/usr/bin/env python3
"""
Shared per-family relationship-role labels for the long-read visualization plots.

`assign_labels(ped)` returns {family: {sample: role}} where role is a short, human-readable
relationship code anchored on the proband(s). The IGV track labels (makeigvpesr.py), the depth
legend (plot_longread_depth.py) and the pedigree nodes (draw_pedigree.py) all import this one
function, so an individual carries the SAME role on every surface and the reader can map a
pedigree symbol straight to its IGV track and depth line.

Roles (numbering is in ped order; the numeric suffix is dropped when a role has one member):
    Pro / Pro1, Pro2...   affected children -- affected individuals that have a parent in the
                          family (falling back to any affected individual if none has a parent)
    Mo / Fa               the proband(s)' mother / father (role wins even if the parent is
                          itself affected; unknown-sex parent -> Par)
    Sib / Sib1, Sib2...   unaffected children sharing >=1 parent with a proband
    MGM MGF PGM PGF       grandparents: maternal/paternal grand-mother/father (unknown sex -> ..P)
    Rel1, Rel2...         anyone else in the family

Ped columns (first six, whitespace- or tab-separated): FamilyID, IndividualID, FatherID,
MotherID, Sex (1=male, 2=female, else unknown), Affected (2=affected). Missing parents are "0"
(also treats "", ".", "-9" as missing). Header lines are skipped. Roles are filename-safe (no
"/"), so they drop straight into IGV symlink names.
"""

import sys
from collections import OrderedDict

_MISSING = {"0", "", ".", "-9"}


def _is_missing(x):
    return x is None or x in _MISSING


def read_ped_families(ped_path):
    """{family: [member, ...]} preserving ped order; member is a dict with keys
    iid, fa, mo, sex, aff (raw strings), idx (0-based order within the family)."""
    fams = OrderedDict()
    with open(ped_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "FamilyID", "family_id")):
                continue
            f = line.split()
            if len(f) < 6:
                continue
            fam, iid, fa, mo, sex, aff = f[:6]
            members = fams.setdefault(fam, [])
            if any(m["iid"] == iid for m in members):     # ped may repeat a row; keep the first
                continue
            members.append({"iid": iid, "fa": fa, "mo": mo, "sex": sex,
                            "aff": aff, "idx": len(members)})
    return fams


def _emit(roles, members, code_of):
    """Assign each member a role from code_of(member), appending a 1-based numeric suffix only
    when a code is shared by more than one member (ped order). Mutates `roles` in place."""
    buckets = OrderedDict()
    for m in sorted(members, key=lambda x: x["idx"]):
        buckets.setdefault(code_of(m), []).append(m)
    for code, ms in buckets.items():
        if len(ms) == 1:
            roles[ms[0]["iid"]] = code
        else:
            for i, m in enumerate(ms, 1):
                roles[m["iid"]] = "{}{}".format(code, i)


def _roles_for_family(members):
    by_iid = {m["iid"]: m for m in members}

    def parent_ids(m):
        ps = []
        for pid in (m["fa"], m["mo"]):
            if not _is_missing(pid) and pid in by_iid:
                ps.append(pid)
        return ps

    roles = {}
    assigned = set()

    affected = [m for m in members if m["aff"] == "2"]
    # probands: affected with a parent in the family; fall back to any affected (lone index)
    candidates = [m for m in affected if parent_ids(m)] or list(affected)
    # an affected individual who is also a parent of another candidate is the parent, not a
    # proband (e.g. an affected mother in a dominant 3-gen family) -> exclude from probands
    parents_of_candidates = set()
    for m in candidates:
        parents_of_candidates.update(parent_ids(m))
    probands = [m for m in candidates if m["iid"] not in parents_of_candidates] or candidates
    _emit(roles, probands, lambda m: "Pro")
    assigned.update(m["iid"] for m in probands)

    # parents of the probands -> Mo / Fa (by sex)
    proband_parent_ids = []
    seen = set()
    for m in probands:
        for pid in parent_ids(m):
            if pid not in seen:
                seen.add(pid)
                proband_parent_ids.append(pid)
    mothers = [by_iid[p] for p in proband_parent_ids
               if by_iid[p]["sex"] == "2" and p not in assigned]
    fathers = [by_iid[p] for p in proband_parent_ids
               if by_iid[p]["sex"] == "1" and p not in assigned]
    parents = [by_iid[p] for p in proband_parent_ids if p not in assigned]
    _emit(roles, parents,
          lambda m: "Mo" if m["sex"] == "2" else ("Fa" if m["sex"] == "1" else "Par"))
    assigned.update(m["iid"] for m in parents)

    # grandparents: parents (in family) of the proband's mother(s) / father(s)
    def grandparents(parent_group):
        out, gseen = [], set()
        for par in parent_group:
            for gid in parent_ids(par):
                if gid not in assigned and gid not in gseen:
                    gseen.add(gid)
                    out.append(by_iid[gid])
        return out

    def gp_code(prefix):
        return lambda m: prefix + ("M" if m["sex"] == "2"
                                   else ("F" if m["sex"] == "1" else "P"))

    maternal_gp = grandparents(mothers)
    paternal_gp = grandparents(fathers)
    _emit(roles, maternal_gp, gp_code("MG"))
    assigned.update(m["iid"] for m in maternal_gp)
    _emit(roles, paternal_gp, gp_code("PG"))
    assigned.update(m["iid"] for m in paternal_gp)

    # unaffected individuals sharing a parent with a proband -> siblings
    proband_parents = set()
    for m in probands:
        proband_parents.update(parent_ids(m))
    sibs = [m for m in members
            if m["iid"] not in assigned and m["aff"] != "2"
            and set(parent_ids(m)) & proband_parents]
    _emit(roles, sibs, lambda m: "Sib")
    assigned.update(m["iid"] for m in sibs)

    # anything left over
    rest = [m for m in members if m["iid"] not in assigned]
    _emit(roles, rest, lambda m: "Rel")
    return roles


def assign_labels(ped_path):
    """{family: {sample: role}} for every family in the ped."""
    return {fam: _roles_for_family(members)
            for fam, members in read_ped_families(ped_path).items()}


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: sample_labels.py <ped>")
    for fam, roles in assign_labels(sys.argv[1]).items():
        for iid, role in sorted(roles.items(), key=lambda kv: kv[1]):
            print("{}\t{}\t{}".format(fam, role, iid))


if __name__ == "__main__":
    main()
