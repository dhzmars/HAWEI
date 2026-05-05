#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

TARGET_RE = re.compile(r"^(L\d+_S\d+_E\d+)__(H\d{3,})$")


@dataclass
class SampleEntry:
    sample: str
    r1: Path
    r2: Path


def parse_args():
    ap = argparse.ArgumentParser(description="Run HAWEI main panels: salmon index -> quant -> CwS")
    ap.add_argument("-o", "--outdir", required=True, help="route directory")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--fq", help="directory containing sample_1.fastq.gz and sample_2.fastq.gz")
    g.add_argument("--fq-manifest", help="TSV with columns: sample, r1, r2")
    ap.add_argument("--samples", default=None, help="optional sample list, comma/space separated")
    ap.add_argument("--panels", default=None, help="optional panel list, comma/space separated; default auto-discover L*_main.fasta")
    ap.add_argument("--nt", default="auto", help="threads; auto -> cpu-2")
    ap.add_argument("--salmon", default="salmon", help="salmon executable [salmon]")
    ap.add_argument("--kmer", type=int, default=31, help="salmon index -k [31]")
    ap.add_argument("--force-index", action="store_true")
    ap.add_argument("--force-quant", action="store_true")
    ap.add_argument("--force-cws", action="store_true")
    return ap.parse_args()


def resolve_nt(s: str) -> int:
    if str(s).lower() == "auto":
        return max(1, (os.cpu_count() or 1) - 2)
    n = int(s)
    if n < 1:
        raise ValueError("--nt must be >= 1")
    return n


def parse_listish(s: Optional[str]) -> Optional[List[str]]:
    if s is None:
        return None
    vals = [x for x in re.split(r"[ ,]+", s.strip()) if x]
    return vals or None


def ensure_dirs(outdir: Path):
    for p in [outdir / "calib" / "idx", outdir / "calib" / "quant", outdir / "calib" / "cws", outdir / "logs"]:
        p.mkdir(parents=True, exist_ok=True)


def parse_manifest(path: Path, wanted: Optional[set[str]]) -> List[SampleEntry]:
    rows: List[SampleEntry] = []
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter='\t')
        need = {"sample", "r1", "r2"}
        if not r.fieldnames or not need.issubset(set(r.fieldnames)):
            raise ValueError(f"Manifest must contain columns: {sorted(need)}")
        for row in r:
            s = row["sample"]
            if wanted and s not in wanted:
                continue
            rows.append(SampleEntry(s, Path(row["r1"]).resolve(), Path(row["r2"]).resolve()))
    if not rows:
        raise ValueError("No samples found in manifest after filtering")
    return rows


def scan_fastq_dir(path: Path, wanted: Optional[set[str]]) -> List[SampleEntry]:
    rows: List[SampleEntry] = []
    seen = set()
    p1s = sorted(path.glob("*_1.fastq.gz")) + sorted(path.glob("*_1.fq.gz"))
    for r1 in p1s:
        sample = re.sub(r"_1\.(fastq|fq)\.gz$", "", r1.name)
        if wanted and sample not in wanted:
            continue
        r2 = None
        for suffix in ["_2.fastq.gz", "_2.fq.gz"]:
            cand = path / f"{sample}{suffix}"
            if cand.exists():
                r2 = cand
                break
        if r2 is not None:
            rows.append(SampleEntry(sample, r1.resolve(), r2.resolve()))
            seen.add(sample)
    if wanted:
        missing = wanted - seen
        if missing:
            raise FileNotFoundError(f"Missing FASTQ pairs for samples: {','.join(sorted(missing))}")
    if not rows:
        raise ValueError(f"No FASTQ pairs found under: {path}")
    return rows


def discover_panels(outdir: Path, wanted: Optional[set[str]]) -> List[str]:
    fasta_dir = outdir / "calib" / "fasta"
    found = sorted(p.stem for p in fasta_dir.glob("L*_main.fasta"))
    if wanted:
        found = [x for x in found if x in wanted]
    if not found:
        raise FileNotFoundError(f"No panel FASTA found under: {fasta_dir}")
    return found


def count_headers(path: Path) -> int:
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


def file_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def index_complete(index_dir: Path) -> bool:
    return file_nonempty(index_dir / "versionInfo.json")


def quant_complete(quant_dir: Path) -> bool:
    if not file_nonempty(quant_dir / "quant.sf"):
        return False
    if not file_nonempty(quant_dir / "aux_info" / "meta_info.json"):
        return False
    if not list((quant_dir / "aux_info").glob("eq_classes.txt*")):
        return False
    return True


def cws_complete(path: Path) -> bool:
    if not file_nonempty(path):
        return False
    with open(path) as f:
        n = 0
        for _ in f:
            n += 1
            if n > 1:
                return True
    return False


def calc_parallel(nt: int, n_tasks: int) -> Tuple[int, int]:
    if n_tasks <= 1:
        return 1, nt
    if nt <= 4:
        jobs = 1
    elif nt <= 12:
        jobs = 2
    elif nt <= 20:
        jobs = 3
    else:
        jobs = 4
    jobs = max(1, min(jobs, n_tasks))
    threads = max(1, nt // jobs)
    return jobs, threads


def open_maybe_gz(path: Path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')


def eq_to_cws(eq_path: Path, out_path: Path):
    with open_maybe_gz(eq_path) as f:
        lines = [x.rstrip('\n') for x in f if x.strip() != ""]
    if len(lines) < 3:
        with open(out_path, 'w') as out:
            out.write("window\thapset\tcount\n")
        return
    n_targets = int(lines[0])
    n_eq = int(lines[1])
    targets = lines[2:2+n_targets]
    eq_lines = lines[2+n_targets:2+n_targets+n_eq]

    parsed: List[Optional[Tuple[str, str]]] = []
    for t in targets:
        m = TARGET_RE.match(t)
        parsed.append((m.group(1), m.group(2)) if m else None)

    agg: Dict[Tuple[str, str], float] = {}
    for line in eq_lines:
        toks = line.split('\t')
        if len(toks) == 1:
            toks = line.split()
        if len(toks) < 2:
            continue
        k = int(toks[0])
        ids = list(map(int, toks[1:1+k]))
        cnt = float(toks[-1])

        windows = []
        haps = []
        ok = True
        for i in ids:
            if i < 0 or i >= len(parsed) or parsed[i] is None:
                ok = False
                break
            w, h = parsed[i]
            windows.append(w)
            haps.append(h)
        if not ok:
            continue
        window_set = "+".join(sorted(set(windows)))
        hapset = ",".join(sorted(set(haps)))
        agg[(window_set, hapset)] = agg.get((window_set, hapset), 0.0) + cnt

    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(["window", "hapset", "count"])
        for (win, hs), c in sorted(agg.items(), key=lambda kv: (-kv[1], kv[0][0], kv[0][1])):
            w.writerow([win, hs, f"{c:.6f}"])


def run_cmd(cmd: List[str], label: str):
    print(f"[RUN] {label}: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main():
    args = parse_args()
    nt = resolve_nt(args.nt)
    print(f"Detected CPU cores: {os.cpu_count() or 1}")
    print(f"--nt {args.nt} resolved to: {nt}")

    outdir = Path(args.outdir).resolve()
    ensure_dirs(outdir)
    sample_filter = set(parse_listish(args.samples) or []) or None
    panel_filter = set(parse_listish(args.panels) or []) or None

    if args.fq_manifest:
        samples = parse_manifest(Path(args.fq_manifest).resolve(), sample_filter)
    else:
        samples = scan_fastq_dir(Path(args.fq).resolve(), sample_filter)

    panels = discover_panels(outdir, panel_filter)
    print(f"Samples: {','.join(s.sample for s in samples)}")
    print(f"Panels: {','.join(panels)}")

    summary_rows = []

    # index
    for panel in panels:
        fa = outdir / "calib" / "fasta" / f"{panel}.fasta"
        idx = outdir / "calib" / "idx" / panel
        idx.mkdir(parents=True, exist_ok=True)
        if count_headers(fa) == 0:
            raise ValueError(f"Empty panel FASTA: {fa}")
        if args.force_index or not index_complete(idx):
            run_cmd([args.salmon, "index", "-t", str(fa), "-i", str(idx), "-k", str(args.kmer)], f"salmon index {panel}")
            idx_status = "BUILT"
        else:
            print(f"[SKIP] index complete: {idx}")
            idx_status = "SKIP"
        summary_rows.append({"panel": panel, "sample": "*index*", "index_status": idx_status, "quant_status": "", "cws_status": "", "index_dir": str(idx), "quant_dir": "", "cws_tsv": ""})

    # quant + cws
    tasks = [(panel, s) for panel in panels for s in samples]
    jobs, threads_per_job = calc_parallel(nt, len(tasks))
    print(f"Parallel jobs: {jobs}; threads/job: {threads_per_job}")

    def worker(panel: str, samp: SampleEntry):
        idx = outdir / "calib" / "idx" / panel
        qdir = outdir / "calib" / "quant" / panel / samp.sample
        qdir.mkdir(parents=True, exist_ok=True)
        cws = outdir / "calib" / "cws" / panel / f"{samp.sample}.CwS.tsv"
        cws.parent.mkdir(parents=True, exist_ok=True)

        quant_status = "SKIP"
        cws_status = "SKIP"

        if args.force_quant or not quant_complete(qdir):
            cmd = [
                args.salmon, "quant",
                "-i", str(idx),
                "-l", "A",
                "-1", str(samp.r1),
                "-2", str(samp.r2),
                "-p", str(threads_per_job),
                "--validateMappings",
                "--dumpEq",
                "--minAssignedFrags", "1",
                "-o", str(qdir),
            ]
            run_cmd(cmd, f"salmon quant {panel} {samp.sample}")
            quant_status = "BUILT"
        else:
            print(f"[SKIP] quant complete: {qdir}")

        eqs = sorted((qdir / "aux_info").glob("eq_classes.txt*"))
        if not eqs:
            raise FileNotFoundError(f"Missing eq_classes under: {qdir / 'aux_info'}")
        eq_path = eqs[0]

        if args.force_cws or not cws_complete(cws):
            eq_to_cws(eq_path, cws)
            if not cws_complete(cws):
                raise ValueError(f"Generated empty CwS: {cws}")
            cws_status = "BUILT"
        else:
            print(f"[SKIP] CwS complete: {cws}")

        return {
            "panel": panel,
            "sample": samp.sample,
            "index_status": "",
            "quant_status": quant_status,
            "cws_status": cws_status,
            "index_dir": str(idx),
            "quant_dir": str(qdir),
            "cws_tsv": str(cws),
        }

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = [ex.submit(worker, panel, samp) for panel, samp in tasks]
        for fu in as_completed(futs):
            summary_rows.append(fu.result())

    summary_path = outdir / "calib" / "run_main_panels.summary.tsv"
    with open(summary_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=["panel", "sample", "index_status", "quant_status", "cws_status", "index_dir", "quant_dir", "cws_tsv"], delimiter='\t')
        w.writeheader()
        for row in sorted(summary_rows, key=lambda x: (x["panel"], x["sample"])):
            w.writerow(row)

    print("Main panel run finished.")
    print(f"Panels discovered: {','.join(panels)}")
    print(f"Index directory: {outdir / 'calib' / 'idx'}")
    print(f"Quant directory: {outdir / 'calib' / 'quant'}")
    print(f"CwS directory: {outdir / 'calib' / 'cws'}")
    print(f"Summary table: {summary_path}")


if __name__ == "__main__":
    main()
