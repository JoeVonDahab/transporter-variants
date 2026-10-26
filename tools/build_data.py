"""Build the static JSON the site reads, from the frozen release files.

Usage:  python -I tools/build_data.py [--release DIR] [--refresh-uniprot]

Writes docs/data/ and copies the release files into docs/downloads/.
The two tracks are written to separate trees and never merged: the native
track carries no third-party predictions (CC BY 4.0), the routed track
inherits CC BY-NC-SA 4.0 from AlphaMissense.
"""
import argparse, gzip, hashlib, json, math, os, re, shutil, sys, time
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
DATA = DOCS / "data"
CACHE = ROOT / "tools" / "cache"
TV = Path(os.path.expanduser("~/projects/TripleVar/outputs/transporters"))

CONS = ["missense", "synonymous", "intronic", "splice_region", "splice_canonical",
        "frameshift", "nonsense", "inframe_indel", "utr3", "utr5", "start_loss",
        "stop_loss", "other"]
EV = ["pathogenic very strong", "pathogenic strong", "pathogenic moderate",
      "pathogenic supporting", "indeterminate", "benign supporting",
      "benign moderate", "benign strong"]
MISSENSE = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")
UNIPROT_FIELDS = ("accession,protein_name,length,sequence,ft_transmem,ft_intramem,"
                  "ft_topo_dom,ft_binding,ft_act_site,ft_site,ft_domain,ft_motif")


def rnd(x, nd):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    return round(float(x), nd)


def sig(x, n=4):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    if x == 0:
        return 0.0
    return float(f"{x:.{n}g}")


def b01(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    if isinstance(x, str):
        return 1 if x.strip().lower() == "true" else 0
    return 1 if bool(x) else 0


def fetch_uniprot(accs, refresh):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / "uniprot.json"
    have = {}
    if path.exists() and not refresh:
        have = json.loads(path.read_text())
    need = [a for a in sorted(set(accs)) if a not in have["entries"]] if have else sorted(set(accs))
    if not have:
        have = {"release": None, "fetched": None, "entries": {}}
    s = requests.Session()
    s.headers["User-Agent"] = "transporter-variants-site-build (github.com/JoeVonDahab)"
    for i in range(0, len(need), 100):
        chunk = need[i:i + 100]
        for attempt in range(4):
            r = s.get("https://rest.uniprot.org/uniprotkb/accessions",
                      params={"accessions": ",".join(chunk), "fields": UNIPROT_FIELDS,
                              "format": "json", "size": 500}, timeout=120)
            if r.status_code == 200:
                break
            time.sleep(3 * (attempt + 1))
        r.raise_for_status()
        have["release"] = r.headers.get("X-UniProt-Release", have["release"])
        for e in r.json()["results"]:
            have["entries"][e["primaryAccession"]] = e
        print(f"  uniprot {min(i + 100, len(need))}/{len(need)}", flush=True)
    if need:
        have["fetched"] = time.strftime("%Y-%m-%d")
        path.write_text(json.dumps(have))
    return have


def features(entry):
    out = []
    for f in entry.get("features", []):
        loc = f.get("location", {})
        try:
            a, b = int(loc["start"]["value"]), int(loc["end"]["value"])
        except (KeyError, TypeError, ValueError):
            continue
        t = f["type"]
        desc = f.get("description") or ""
        if t == "Binding site" and f.get("ligand"):
            desc = f["ligand"].get("name", desc)
        out.append([t, a, b, desc])
    return out


def protein_name(entry):
    pd_ = entry.get("proteinDescription", {})
    rec = pd_.get("recommendedName") or (pd_.get("submissionNames") or [{}])[0]
    return rec.get("fullName", {}).get("value", "")


def parse_calibration_v11(md):
    """CALIBRATION_CHECK.md from release 1.1.0 on: score thresholds plus one validation table per track."""
    out = {"method": "pejaver2022", "prior": 0.0441, "thresholds": [], "tables": {}}
    sections = re.split(r"^#{2,3} ", md, flags=re.M)
    for sec in sections:
        head, _, body = sec.partition("\n")
        head = head.strip()
        if head.startswith("Score thresholds"):
            for line in body.splitlines():
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                if len(cells) >= 11 and cells[0] not in ("Unit", "") and not cells[0].startswith("---"):
                    out["thresholds"].append({"unit": cells[0], "serves": cells[1], "labelled": cells[2],
                                              "levels": cells[3:10], "max": cells[10]})
            continue
        m = re.match(r"(routed|native) track, (ClinGen calibration|missense|for comparison)", head)
        if not m:
            continue
        key = f"{m.group(1)}_{'all' if m.group(2).startswith('ClinGen') else 'missense' if m.group(2) == 'missense' else 'old'}"
        rows = []
        for line in body.splitlines():
            r = re.match(r"^\|\s*([a-z ]+?)(?:\s*\(([^)]*)\))?\s*\|\s*([\d,]+)\s*\|\s*([\d.]+)%\s*\|\s*([\d.e+-]+)\s*\[([\d.e+-]+),\s*([\d.e+-]+)\]\s*\|\s*([\d.]+)\s*\|\s*(.*?)\s*\|$", line)
            if r:
                rows.append({"level": r.group(1).strip(), "code": r.group(2), "n": int(r.group(3).replace(",", "")),
                             "observed": float(r.group(4)), "lr": float(r.group(5)), "lo": float(r.group(6)),
                             "hi": float(r.group(7)), "posterior": float(r.group(8)), "flag": r.group(9)})
        out["tables"][key] = rows
    m = re.search(r"([\d,]+) labelled out-of-sample variants", md)
    out["n"] = int(m.group(1).replace(",", "")) if m else None
    return out


def parse_calibration(md):
    rows, base, n_total, n_genes = [], None, None, None
    for line in md.splitlines():
        m = re.match(r"^\|\s*([a-z ]+?)\s*\|\s*([\d,]+)\s*\|\s*([\d.]+)%\s*\|\s*([\d.]+)x\s*\|$", line)
        if m:
            rows.append({"bin": m.group(1), "n": int(m.group(2).replace(",", "")),
                         "observed": float(m.group(3)), "vs_base": float(m.group(4))})
    m = re.search(r"([\d,]+) labelled out-of-sample variants across ([\d,]+) genes; "
                  r"the base rate is ([\d.]+)% pathogenic", md)
    if m:
        n_total, n_genes, base = int(m.group(1).replace(",", "")), int(m.group(2)), float(m.group(3))
    return {"rows": rows, "n": n_total, "genes": n_genes, "base_rate": base}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", default=str(TV / "release"))
    ap.add_argument("--refresh-uniprot", action="store_true")
    a = ap.parse_args()
    REL = Path(a.release)

    print("reading release", flush=True)
    nat = pd.read_csv(REL / "transporter_variants_native.csv.gz", low_memory=False)
    rou = pd.read_csv(REL / "transporter_variants_routed.csv.gz", low_memory=False)
    assert len(nat) == len(rou) == 349056, (len(nat), len(rou))
    assert (nat.variation_id.values == rou.variation_id.values).all()
    assert set(nat.evidence.unique()) <= set(EV) and set(rou.evidence.unique()) <= set(EV)
    assert set(nat.consequence.unique()) <= set(CONS)
    nat["family"] = nat.family.fillna("unassigned")
    rou["family"] = rou.family.fillna("unassigned")

    # Which branch scored each variant. The release drops this column, but it matters for missense:
    # a missense variant whose protein change could not be mapped onto the UniProt sequence was scored
    # by the consequence-class branch, which knows only that it is missense, not where or what.
    # from 1.0.2 on, rescued variants moved to the missense branch; their scores live in the
    # versioned file, while published_scores.csv stays frozen as the 1.0.1 manuscript input
    scores_file = TV / "published_scores_v1.0.2.csv"
    if not scores_file.exists():
        scores_file = TV / "published_scores.csv"
    print(f"branch from {scores_file.name}")
    br = pd.read_csv(scores_file, usecols=["variation_id", "branch"])
    br = dict(zip(br.variation_id, br.branch))
    for df in (nat, rou):
        df["class_level"] = ((df.consequence == "missense") &
                             (df.variation_id.map(br) == "consequence")).astype(int)
    print(f"missense scored by the consequence branch: {int(nat.class_level.sum()):,}")

    # gene -> UniProt accession, exactly as 114_add_layers.py mapped them
    prot = pd.concat([pd.read_csv(TV / "proteins_dms.csv"),
                      pd.read_csv(TV / "proteins_all.csv").dropna(subset=["gene"])])
    prot = prot.drop_duplicates("gene")
    g2a = dict(zip(prot.gene.astype(str), prot.uniprot_id))
    g2seq = dict(zip(prot.gene.astype(str), prot.sequence))
    genes = sorted(nat.gene.astype(str).unique())
    missing_acc = [g for g in genes if g not in g2a]
    print(f"{len(genes)} genes, {len(missing_acc)} without an accession: {missing_acc[:20]}")

    up = fetch_uniprot([g2a[g] for g in genes if g in g2a], a.refresh_uniprot)

    for sub in ("native", "routed", "protein", "lookup/rs", "lookup/vid"):
        p = DATA / sub
        if p.exists():
            shutil.rmtree(p)
        p.mkdir(parents=True)

    index, report = [], []
    rs_lookup, vid_lookup = defaultdict(dict), defaultdict(dict)
    # repeated strings are stored once, in index.json, and referenced by position
    CLINSIG = sorted(nat.clinsig.dropna().astype(str).unique())
    REVIEW = sorted(nat.review_status.dropna().astype(str).unique())
    cs_i = {s: i for i, s in enumerate(CLINSIG)}
    rv_i = {s: i for i, s in enumerate(REVIEW)}

    def split_hgvs(h, prefixes):
        # "NM_006996.3(SLC19A2):c.725del (p.Pro242fs)" -> prefix index 0, "c.725del (p.Pro242fs)"
        h = str(h)
        k = h.find(":")
        if k < 0:
            return [-1, h]
        p = h[:k + 1]
        if p not in prefixes:
            prefixes[p] = len(prefixes)
        return [prefixes[p], h[k + 1:]]
    groups_n = dict(tuple(nat.groupby("gene", sort=False)))
    groups_r = dict(tuple(rou.groupby("gene", sort=False)))

    for g in genes:
        dn, dr = groups_n[g], groups_r[g]
        fam = dn.family.iloc[0]
        acc = g2a.get(g)
        seq = g2seq.get(g) if acc else None
        e = up["entries"].get(acc) if acc else None

        # ---- protein meta: track-independent, no predictions
        ref_ok = ref_bad = 0
        pos_fields = []
        for pv in dn.protein_variant.astype(str):
            m = MISSENSE.match(pv)
            if not m:
                pos_fields.append(None)
                continue
            pos = int(m.group(2))
            pos_fields.append(pos)
            if seq and 1 <= pos <= len(seq) and seq[pos - 1] == m.group(1):
                ref_ok += 1
            else:
                ref_bad += 1
        up_seq = e.get("sequence", {}).get("value") if e else None
        meta = {
            "gene": g, "family": fam, "uniprot": acc,
            "name": protein_name(e) if e else "",
            "length": len(seq) if seq else None,
            "sequence": seq,
            "uniprot_sequence_matches": (up_seq == seq) if (up_seq and seq) else None,
            "features": features(e) if e and up_seq == seq else [],
            "ref_match": [ref_ok, ref_bad],
        }
        (DATA / "protein" / f"{g}.json").write_text(json.dumps(meta, separators=(",", ":")))

        # ---- shared identity columns, then each track's predictions.
        # hgvs is stored as [prefix index, remainder]; clinsig and review status
        # as positions in index.json's dictionaries.
        def base(t, pos, prefixes):
            rs = int(t.rsid) if pd.notna(t.rsid) and int(t.rsid) > 0 else None
            return [int(t.variation_id), rs, *split_hgvs(t.hgvs, prefixes),
                    CONS.index(t.consequence),
                    None if pd.isna(t.protein_variant) else str(t.protein_variant),
                    None if pd.isna(t.clinsig) else cs_i[str(t.clinsig)],
                    None if pd.isna(t.review_status) else rv_i[str(t.review_status)],
                    b01(t.in_training), EV.index(t.evidence), rnd(t.score, 4), sig(t.lr),
                    rnd(t.activity_pred, 3), rnd(t.mechanism_pred, 3),
                    b01(t.layers_in_sample), pos, int(t.class_level),
                    None if pd.isna(getattr(t, "acmg_code", None)) else str(t.acmg_code),
                    None if pd.isna(getattr(t, "evidence_v1", None)) or t.evidence_v1 not in EV else EV.index(t.evidence_v1)]
        COLS = ["variation_id", "rsid", "hgvs_prefix", "hgvs", "consequence", "protein_variant",
                "clinsig", "review_status", "in_training", "evidence", "score", "lr",
                "activity_pred", "mechanism_pred", "layers_in_sample", "pos", "class_level",
                "acmg_code", "evidence_v1"]

        # ---- native track (CC BY 4.0)
        rows, pre = [], {}
        for t, pos in zip(dn.itertuples(index=False), pos_fields):
            rows.append(base(t, pos, pre))
            rs = rows[-1][1]
            if rs:
                lst = rs_lookup[f"{rs % 100:02d}"].setdefault(str(rs), [])
                if g not in lst:
                    lst.append(g)
            vid_lookup[f"{int(t.variation_id) % 100:02d}"][str(int(t.variation_id))] = g
        (DATA / "native" / f"{g}.json").write_text(json.dumps({
            "gene": g, "track": "native", "license": "CC BY 4.0",
            "hgvs_prefixes": list(pre), "cols": COLS, "rows": rows}, separators=(",", ":")))

        # ---- routed track (CC BY-NC-SA 4.0)
        rows_r, pre_r = [], {}
        for t, pos in zip(dr.itertuples(index=False), pos_fields):
            rows_r.append(base(t, pos, pre_r) + [
                1 if t.track == "plus" else 0, b01(t.has_alphamissense),
                rnd(t.v2_native, 4), rnd(t.v2_plus, 4), rnd(t.alphamissense, 4)])
        (DATA / "routed" / f"{g}.json").write_text(json.dumps({
            "gene": g, "track": "routed", "license": "CC BY-NC-SA 4.0",
            "hgvs_prefixes": list(pre_r),
            "cols": COLS + ["served_by_plus", "has_alphamissense", "v2_native", "v2_plus",
                            "alphamissense"],
            "rows": rows_r}, separators=(",", ":")))

        bins_n = Counter(dn.evidence)
        bins_r = Counter(dr.evidence)
        tx = sorted({h.split("(")[0].split(".")[0] for h in dn.hgvs.astype(str)
                     if h.startswith(("NM_", "NR_", "XM_"))})
        index.append([g, fam, acc, meta["name"], len(dn), int((dn.consequence == "missense").sum()),
                      [bins_n.get(b, 0) for b in EV], [bins_r.get(b, 0) for b in EV],
                      int(dn.in_training.sum()), tx])
        report.append((g, ref_ok, ref_bad, meta["uniprot_sequence_matches"], len(meta["features"])))

    for kind, lk in (("rs", rs_lookup), ("vid", vid_lookup)):
        for k, v in lk.items():
            (DATA / "lookup" / kind / f"{k}.json").write_text(json.dumps(v, separators=(",", ":")))

    # ---- landing-page numbers, all computed from the files themselves
    by_cons = []
    for c in CONS:
        s = nat[nat.consequence == c]
        by_cons.append([c, len(s), int(s.has_alphamissense.sum())])
    fam_counts = nat.groupby("family").agg(genes=("gene", "nunique"), variants=("gene", "size"))
    cal_md = (REL / "CALIBRATION_CHECK.md").read_text()
    cal = parse_calibration_v11(cal_md) if "Pejaver" in cal_md else parse_calibration(cal_md)
    if cal.get("method") == "pejaver2022":
        for k in ("routed_all", "native_all", "routed_missense", "native_missense"):
            assert cal["tables"].get(k), f"calibration table {k} not parsed"
        assert len(cal["thresholds"]) == 3, "threshold table not parsed"
    cff = (REL / "CITATION.cff").read_text()
    version = re.search(r"^version:\s*(\S+)", cff, re.M).group(1)
    released = re.search(r"^date-released:\s*(\S+)", cff, re.M).group(1)
    summary = {
        "version": version, "released": released,
        "variants": len(nat), "genes": len(genes),
        "missense": int((nat.consequence == "missense").sum()),
        "in_training": int(nat.in_training.sum()),
        "has_alphamissense": int(nat.has_alphamissense.sum()),
        "activity": int(nat.activity_pred.notna().sum()),
        "mechanism_native": int(nat.mechanism_pred.notna().sum()),
        "mechanism_routed": int(rou.mechanism_pred.notna().sum()),
        "routed_plus": int((rou.track == "plus").sum()),
        "class_level_missense": int(nat.class_level.sum()),
        "acmg_native": nat.acmg_code.value_counts().to_dict() if "acmg_code" in nat else {},
        "acmg_routed": rou.acmg_code.value_counts().to_dict() if "acmg_code" in rou else {},
        "by_consequence": by_cons,
        "families": [[f, int(r.genes), int(r.variants)] for f, r in
                     fam_counts.sort_values("variants", ascending=False).iterrows()],
        "evidence_native": [int((nat.evidence == b).sum()) for b in EV],
        "evidence_routed": [int((rou.evidence == b).sum()) for b in EV],
        "calibration": cal,
        "uniprot_release": up.get("release"),
    }

    # ---- downloads: each track in its own folder with its own licence
    dl = DOCS / "downloads"
    for sub, files in (("native", ["transporter_variants_native.csv.gz", "LICENSE-native"]),
                       ("routed", ["transporter_variants_routed.csv.gz", "LICENSE-routed"]),
                       ("", ["README.md", "DATA_DICTIONARY.md", "CALIBRATION_CHECK.md",
                             "CITATION.cff", "CHECKSUMS.txt", "CHANGELOG.md", "EVIDENCE_1.1.0.md",
                             "COMPLETION_115.md"])):
        (dl / sub).mkdir(parents=True, exist_ok=True)
        for f in files:
            if (REL / f).exists():
                shutil.copy2(REL / f, dl / sub / f)
    sums = dict(l.split()[::-1] for l in (REL / "CHECKSUMS.txt").read_text().splitlines() if l.strip())
    for sub, f in (("native", "transporter_variants_native.csv.gz"),
                   ("routed", "transporter_variants_routed.csv.gz")):
        got = sha256(dl / sub / f)
        assert got == sums[f], f"checksum mismatch for {f}"
    summary["downloads"] = {
        "native": {"bytes": (dl / "native/transporter_variants_native.csv.gz").stat().st_size,
                   "sha256": sums["transporter_variants_native.csv.gz"]},
        "routed": {"bytes": (dl / "routed/transporter_variants_routed.csv.gz").stat().st_size,
                   "sha256": sums["transporter_variants_routed.csv.gz"]},
    }
    (DATA / "index.json").write_text(json.dumps(
        {"cons": CONS, "ev": EV, "clinsig": CLINSIG, "review": REVIEW, "summary": summary,
         "cols": ["gene", "family", "uniprot", "name", "n", "n_missense", "bins_native",
                  "bins_routed", "n_in_training", "transcripts"],
         "genes": index}, separators=(",", ":")))

    # ---- report
    tot_ok = sum(r[1] for r in report); tot_bad = sum(r[2] for r in report)
    print(f"\nmissense reference residue matches the model's sequence: {tot_ok:,} ok, {tot_bad:,} not")
    worst = sorted(report, key=lambda r: -(r[2] / max(1, r[1] + r[2])))[:12]
    print("worst genes:", [(r[0], r[1], r[2]) for r in worst])
    print("uniprot sequence != model sequence:", [r[0] for r in report if r[3] is False])
    print("no uniprot entry:", [r[0] for r in report if r[3] is None])
    print("unassigned family:", sorted(nat[nat.family == "unassigned"].gene.unique().tolist()))
    print("genes with no topology features:", sum(1 for r in report if r[4] == 0))
    print("rsid <= 0:", int((nat.rsid <= 0).sum()))
    print("mechanism outside [0,1]:", int(((nat.mechanism_pred < 0) | (nat.mechanism_pred > 1)).sum()),
          "of", int(nat.mechanism_pred.notna().sum()))
    sz = sum(p.stat().st_size for p in DATA.rglob("*.json"))
    print(f"data size {sz / 1e6:.1f} MB in {sum(1 for _ in DATA.rglob('*.json'))} files")


if __name__ == "__main__":
    main()
