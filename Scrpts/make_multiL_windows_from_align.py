#!/usr/bin/env python3
import argparse
import csv
import os
import re
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import List, Tuple


def parse_args():
    ap = argparse.ArgumentParser(
        description="Build HAWEI windows and main panels from an aligned full-length hap FASTA"
    )
    ap.add_argument("-i", "--input", required=True, help="aligned full-length hap FASTA")
    ap.add_argument("-o", "--outdir", required=True, help="output route directory")
    ap.add_argument("--len", default="300,400,800,1200", help="window lengths, comma or space separated")
    ap.add_argument("--min-ungapped", type=int, default=80, help="minimum ungapped nt required to keep a hap window [80]")
    ap.add_argument("--nt", default="auto", help="threads placeholder for interface consistency; auto -> cpu-2")
    return ap.parse_args()


def resolve_nt(s: str) -> int:
    if str(s).lower() == "auto":
        return max(1, (os.cpu_count() or 1) - 2)
    n = int(s)
    if n < 1:
        raise ValueError("--nt must be >= 1")
    return n


def parse_lengths(s: str) -> List[int]:
    vals = [x for x in re.split(r"[ ,]+", s.strip()) if x]
    lens = [int(x) for x in vals]
    if not lens:
        raise ValueError("No window lengths provided")
    if any(L < 1 for L in lens):
        raise ValueError("Window lengths must be positive")
    return lens


def read_fasta(path: Path) -> OrderedDict:
    seqs = OrderedDict()
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                # enforce first whitespace-delimited token as hap label
                name = line[1:].split()[0]
                if name in seqs:
                    raise ValueError(f"Duplicate FASTA header token: {name}")
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    if len(seqs) < 2:
        raise ValueError("Need at least two aligned hap sequences")
    return seqs


def ensure_aligned(seqs: OrderedDict) -> int:
    lens = {len(v) for v in seqs.values()}
    if len(lens) != 1:
        raise ValueError(f"All aligned sequences must have identical length, got lengths: {sorted(lens)}")
    return next(iter(lens))


SAFE_HAP_RE = re.compile(r"^H\d{3,}$")


def assign_hap_ids(labels: List[str]) -> OrderedDict:
    mapping = OrderedDict()
    for i, lab in enumerate(labels, start=1):
        hid = f"H{i:03d}"
        if SAFE_HAP_RE.fullmatch(lab):
            # still map original to generated internal IDs to avoid collisions with output parsing
            pass
        mapping[hid] = lab
    return mapping


WINDOW_RE = re.compile(r"^(L\d+_S\d+_E\d+)__(H\d{3,})$")


def format_win_id(L: int, start: int, end: int) -> str:
    return f"L{L}_S{start:04d}_E{end:04d}"


def regular_and_tail_anchor_starts(aln_len: int, L: int) -> List[Tuple[int, int]]:
    starts = []
    s = 1
    while s + L - 1 <= aln_len:
        starts.append((s, s + L - 1))
        s += L
    if aln_len % L != 0:
        ts = aln_len - L + 1
        te = aln_len
        if (ts, te) not in starts:
            starts.append((ts, te))
    return starts


def main():
    args = parse_args()
    nt = resolve_nt(args.nt)
    print(f"Detected CPU cores: {os.cpu_count() or 1}")
    print(f"--nt {args.nt} resolved to: {nt}")

    in_fa = Path(args.input).resolve()
    outdir = Path(args.outdir).resolve()
    build_dir = outdir / "build"
    lists_dir = outdir / "calib" / "lists"
    fasta_dir = outdir / "calib" / "fasta"
    logs_dir = outdir / "logs"
    for p in [build_dir, lists_dir, fasta_dir, logs_dir]:
        p.mkdir(parents=True, exist_ok=True)

    seqs = read_fasta(in_fa)
    aln_len = ensure_aligned(seqs)
    lengths = parse_lengths(args.len)
    hap_map = assign_hap_ids(list(seqs.keys()))
    seq_by_hid = OrderedDict((hid, seqs[lab]) for hid, lab in hap_map.items())

    hap_map_tsv = build_dir / "hap_map.tsv"
    with open(hap_map_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(["hap_id", "hap_label", "aligned_length", "ungapped_length"])
        for hid, lab in hap_map.items():
            seq = seqs[lab]
            w.writerow([hid, lab, len(seq), len(seq.replace('-', ''))])

    all_fa = build_dir / "windows_multiL.clean.fasta"
    panel_handles = {}
    windows_by_len = defaultdict(list)
    panel_seq_counts = defaultdict(int)
    panel_window_counts = defaultdict(int)

    with open(all_fa, "w") as all_out:
        for L in lengths:
            pfa = fasta_dir / f"L{L}_main.fasta"
            panel_handles[L] = open(pfa, "w")
            seen_windows = set()
            for start, end in regular_and_tail_anchor_starts(aln_len, L):
                win_id = format_win_id(L, start, end)
                wrote_any = False
                for hid, seq in seq_by_hid.items():
                    frag = seq[start - 1:end]
                    ungapped = frag.replace('-', '')
                    if len(ungapped) < args.min_ungapped:
                        continue
                    header = f"{win_id}__{hid}"
                    all_out.write(f">{header}\n")
                    panel_handles[L].write(f">{header}\n")
                    for i in range(0, len(ungapped), 80):
                        chunk = ungapped[i:i+80]
                        all_out.write(chunk + "\n")
                        panel_handles[L].write(chunk + "\n")
                    wrote_any = True
                    panel_seq_counts[L] += 1
                if wrote_any and win_id not in seen_windows:
                    seen_windows.add(win_id)
                    windows_by_len[L].append(win_id)
                    panel_window_counts[L] += 1
        for h in panel_handles.values():
            h.close()

    for L in lengths:
        list_path = lists_dir / f"L{L}_main.list"
        with open(list_path, "w") as f:
            for win in windows_by_len[L]:
                f.write(win + "\n")

    summary_path = build_dir / "build_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(["panel", "n_windows", "n_sequences", "alignment_length", "n_haps", "input_fasta", "all_windows_fasta", "panel_fasta", "panel_list"])
        for L in lengths:
            w.writerow([
                f"L{L}_main",
                panel_window_counts[L],
                panel_seq_counts[L],
                aln_len,
                len(hap_map),
                str(in_fa),
                str(all_fa),
                str(fasta_dir / f"L{L}_main.fasta"),
                str(lists_dir / f"L{L}_main.list"),
            ])

    print("Build finished.")
    print(f"Input alignment: {in_fa}")
    print(f"Hap map: {hap_map_tsv}")
    print(f"All windows FASTA: {all_fa}")
    print(f"Panel FASTA directory: {fasta_dir}")
    print(f"Panel list directory: {lists_dir}")
    print(f"Build summary: {summary_path}")


if __name__ == "__main__":
    main()
