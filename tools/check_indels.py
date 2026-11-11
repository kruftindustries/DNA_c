#!/usr/bin/env python3
"""Inventory of indel-keyed entries across the data the report is built from.

Array exports write an indel genotype as D/I letters (DD, DI, II); ClinPGx
spells the same genotypes with allele sequences and "del"; ClinVar holds
indels in VCF spelling (anchor base kept). This lists, for a data directory:

  * ClinPGx rsID annotations whose genotype rows are not plain base pairs,
    with the spelling pattern, evidence level and whether the engine can
    derive D/I keys for them (exactly one `del` form and one sequence form);
  * ClinVar indel records at a reportable significance, by rsID, flagging
    rsIDs that carry more than one indel record (which the engine declines
    to call from a D/I genotype);
  * curated-table keys that can never match: non-base keys that are not
    D/I spellings, or base keys at an rsID ClinPGx treats as an indel.

    tools/check_indels.py [--data DIR] [--tables c/src]

Re-run after every data update; the engine's rules live in
c/src/gh_pharmgkb.c (derive_indel_keys) and c/src/gh_clinvar.c.
"""
import argparse
import collections
import csv
import pathlib
import re
import sys

csv.field_size_limit(1 << 30)
SIMPLE = re.compile(r"^[ACGT]{1,2}$")
SEQ = re.compile(r"^[ACGT]+$")
REPORTABLE_PGX = {"1A", "1B", "2A", "2B"}


def derivable(rows):
    """Mirror of derive_indel_keys: every row must map to D/I letters."""
    has_del = any("del" in r.lower() for r in rows)
    if not has_del:
        return False
    for r in rows:
        alleles = r.split("/") if "/" in r else ([r[: len(r) // 2], r[len(r) // 2:]] if len(r) % 2 == 0 and r else [r])
        for a in alleles:
            if a.lower() == "del" or SEQ.match(a):
                continue
            return False
    return True


def clinpgx(data):
    ann_path = data / "clinical_annotations.tsv"
    all_path = data / "clinical_ann_alleles.tsv"
    if not ann_path.exists() or not all_path.exists():
        print("ClinPGx tables not present; run gh-data pharmgkb")
        return {}
    ann = {}
    with ann_path.open(newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            v = r["Variant/Haplotypes"]
            if v.startswith("rs") and "," not in v:
                ann[r.get("Clinical Annotation ID") or r.get("Summary Annotation ID")] = (
                    v, r["Gene"], r["Level of Evidence"], r["Drug(s)"])
    rows = collections.defaultdict(set)
    with all_path.open(newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            a = ann.get(r.get("Clinical Annotation ID") or r.get("Summary Annotation ID"))
            if a:
                rows[a].add(r["Genotype/Allele"])
    indel = {a: g for a, g in rows.items() if any(not SIMPLE.match(x) for x in g)}
    print(f"ClinPGx: {len(ann)} rsID annotations, {len(indel)} with indel-form genotype rows")
    by_rsid = {}
    for (rsid, gene, level, drugs), g in sorted(indel.items(), key=lambda t: (t[0][2], t[0][1])):
        ok = derivable(g)
        by_rsid.setdefault(rsid, ok)
        flag = "D/I keys derived" if ok else "NOT matchable from D/I"
        mark = "*" if level in REPORTABLE_PGX else " "
        shown = [a if len(a) <= 14 else a[:6] + "…" + a[-6:] for a in sorted(g)[:3]]
        print(f"  {mark} {rsid:<12} {gene:<8} {level:<3} {flag:<24} {shown}  {drugs[:40]}")
    print("  (* = evidence level the report shows)")
    return by_rsid


def clinvar(data):
    path = data / "clinvar_alleles.tsv"
    if not path.exists():
        print("ClinVar extract not present; run gh-data clinvar")
        return
    per_rsid = collections.defaultdict(set)
    with_rsid = False
    n = no_alleles = no_rsid = mnv = 0
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        with_rsid = "rsid" in (reader.fieldnames or [])
        for r in reader:
            ref, alt = r["ref"], r["alt"]
            if len(ref) == 1 and len(alt) == 1:
                continue
            sig = r["clinical_significance"].lower()
            if not any(k in sig for k in ("pathogenic", "risk factor", "drug response", "protective")):
                continue
            # variant_summary writes "na" where it has no VCF alleles (large
            # or complex changes); the engine cannot read those from D/I.
            if not (SEQ.match(ref) and SEQ.match(alt)):
                no_alleles += 1
                continue
            if len(ref) == len(alt):
                mnv += 1
                continue
            n += 1
            rsid = r.get("rsid") or ""
            if not rsid:
                no_rsid += 1
                continue
            per_rsid[rsid].add((ref, alt))
    ambiguous = {k: v for k, v in per_rsid.items() if len(v) > 1}
    print(f"ClinVar: {n} deletion/insertion records at a reported significance "
          f"({no_alleles} more have no VCF alleles, {mnv} are MNVs: neither is readable from D/I)")
    if not with_rsid:
        print("    no rsid column: re-run gh-data clinvar for indel matching by rsID")
    else:
        print(f"    {len(per_rsid)} rsIDs, {no_rsid} records without one (matched by position only); "
              f"{len(ambiguous)} rsIDs carry more than one distinct indel and are not called from D/I")
        for k, v in list(sorted(ambiguous.items(), key=lambda t: int(t[0][2:])))[:8]:
            v = sorted(v)
            print(f"    {k}: {len(v)} alleles, e.g. {v[0][0]}>{v[0][1]} and {v[1][0]}>{v[1][1]}")


def tables(src, pgx_indels):
    """Curated keys in the generated tables that cannot be met by any file."""
    text = ""
    for p in sorted(src.glob("gh_*_data.c")):
        text += p.read_text(errors="replace")
    # gh_snp_variant rows: {"GENOTYPE", "status", "description", magnitude}
    keys = collections.Counter(re.findall(r'\{"([^"]+)",\s*"[^"]*",\s*"[^"]*",\s*\d+\}', text))
    odd = {k: c for k, c in keys.items() if not SIMPLE.match(k) and not re.fullmatch(r"[DI]{1,2}", k)
           and k not in ("Unknown",) and not k.startswith("*") and not k.startswith("rs")}
    print(f"Curated tables: {len(odd)} genotype-like keys that are neither bases nor D/I: "
          + ", ".join(f"{k}({c})" for k, c in sorted(odd.items())[:12]))
    for rsid, ok in sorted(pgx_indels.items()):
        m = re.search(r"\.rsid = \"%s\"" % rsid, text)
        if m:
            print(f"  curated entry at ClinPGx indel {rsid}: check its keys are D/I spellings "
                  f"({'derivable' if ok else 'ClinPGx itself not derivable'})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--tables", default="c/src")
    a = ap.parse_args()
    data = pathlib.Path(a.data)
    pgx = clinpgx(data)
    clinvar(data)
    tables(pathlib.Path(a.tables), pgx)
    return 0


if __name__ == "__main__":
    sys.exit(main())
