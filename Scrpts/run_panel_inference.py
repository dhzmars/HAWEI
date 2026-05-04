#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import os
import random
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

TARGET_RE = re.compile(r"^(L\d+_S\d+_E\d+)__(H\d{3,})$")


def parse_args():
    ap = argparse.ArgumentParser(description="Run HAWEI panel inference: ALL3 merge -> EM -> panel selection -> EC bootstrap")
    ap.add_argument("-o", "--outdir", required=True, help="route directory")
    ap.add_argument("--samples", default=None, help="optional sample list, comma/space separated; default auto-detect from CwS")
    ap.add_argument("--panels", default=None, help="optional panel list, comma/space separated; default auto-detect from CwS")
    ap.add_argument("--nt", default="auto", help="threads; auto -> cpu-2")
    ap.add_argument("--em-boot", type=int, default=2000, help="CwS-level bootstrap replicates for panel screening [2000]")
    ap.add_argument("--ec-boot", type=int, default=2000, help="eq_classes-level bootstrap replicates for final panel [2000]")
    ap.add_argument("--seed", type=int, default=1, help="random seed [1]")
    ap.add_argument("--stop-after", choices=["select", "ecboot"], default="ecboot")
    ap.add_argument("--panel", default=None, help="manual panel override for EC bootstrap")
    ap.add_argument("--ecboot-top", type=int, default=1, help="if no --panel, run EC bootstrap for top N panels [1]")
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
    for p in [outdir / "calib" / "em", outdir / "calib" / "ecboot", outdir / "calib" / "summary", outdir / "logs"]:
        p.mkdir(parents=True, exist_ok=True)


def read_hap_map(path: Path) -> List[Tuple[str, str]]:
    rows = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            rows.append((row['hap_id'], row['hap_label']))
    if not rows:
        raise ValueError(f"Empty hap map: {path}")
    return rows


def discover_panels(outdir: Path, wanted: Optional[set[str]]) -> List[str]:
    cws_root = outdir / "calib" / "cws"
    found = sorted(p.name for p in cws_root.iterdir() if p.is_dir()) if cws_root.exists() else []
    if wanted:
        found = [x for x in found if x in wanted]
    if not found:
        raise FileNotFoundError(f"No panel CwS directories found under: {cws_root}")
    return found


def discover_samples(outdir: Path, panel: str, wanted: Optional[set[str]]) -> List[str]:
    panel_dir = outdir / "calib" / "cws" / panel
    samples = []
    for p in sorted(panel_dir.glob("*.CwS.tsv")):
        if p.name == "ALL3.CwS.tsv":
            continue
        s = p.stem.replace('.CwS', '')
        if wanted and s not in wanted:
            continue
        samples.append(s)
    if wanted:
        missing = wanted - set(samples)
        if missing:
            raise FileNotFoundError(f"Missing CwS files for samples under {panel_dir}: {','.join(sorted(missing))}")
    if not samples:
        raise FileNotFoundError(f"No sample CwS found under: {panel_dir}")
    return samples


def read_cws_rows(path: Path) -> List[Tuple[str, Tuple[str, ...], int]]:
    rows = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            hs = tuple(sorted(x for x in row['hapset'].split(',') if x))
            rows.append((row['window'], hs, int(round(float(row['count'])))))
    return rows


def write_cws_rows(path: Path, rows: List[Tuple[str, Tuple[str, ...], int]]):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['window', 'hapset', 'count'])
        for win, hs, c in sorted(rows, key=lambda x: (-x[2], x[0], x[1])):
            w.writerow([win, ','.join(hs), c])


def merge_all3(panel_dir: Path, samples: List[str]) -> Path:
    agg: Dict[Tuple[str, Tuple[str, ...]], int] = {}
    for s in samples:
        for win, hs, c in read_cws_rows(panel_dir / f"{s}.CwS.tsv"):
            agg[(win, hs)] = agg.get((win, hs), 0) + c
    rows = [(win, hs, c) for (win, hs), c in agg.items()]
    out = panel_dir / 'ALL3.CwS.tsv'
    write_cws_rows(out, rows)
    return out


def multinomial_resample(rows: List[Tuple[str, Tuple[str, ...], int]], rng: random.Random) -> List[Tuple[str, Tuple[str, ...], int]]:
    counts = [c for _, _, c in rows]
    total = sum(counts)
    if total <= 0:
        return [(w, hs, 0) for (w, hs, _) in rows]
    probs = [c / total for c in counts]
    cum = []
    s = 0.0
    for p in probs:
        s += p
        cum.append(s)
    new_counts = [0] * len(rows)
    for _ in range(total):
        x = rng.random()
        lo, hi = 0, len(cum) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if x <= cum[mid]:
                hi = mid
            else:
                lo = mid + 1
        new_counts[lo] += 1
    return [(rows[i][0], rows[i][1], new_counts[i]) for i in range(len(rows))]


def em_pi_from_rows(rows: List[Tuple[str, Tuple[str, ...], int]], haps: List[str], max_iter: int = 5000, tol: float = 1e-12, eps: float = 1e-15) -> Dict[str, float]:
    if not rows:
        raise ValueError("No usable rows for EM")
    pi = {h: 1.0 / len(haps) for h in haps}
    haplists = [list(hs) for _, hs, _ in rows]
    counts = [c for _, _, c in rows]
    for _ in range(max_iter):
        exp = {h: 0.0 for h in haps}
        for hs_list, c in zip(haplists, counts):
            denom = sum(pi[h] for h in hs_list)
            if denom <= 0:
                share = c / len(hs_list)
                for h in hs_list:
                    exp[h] += share
            else:
                for h in hs_list:
                    exp[h] += c * (pi[h] / denom)
        total = sum(exp.values())
        if total <= 0:
            pi_new = {h: 1.0 / len(haps) for h in haps}
        else:
            pi_new = {h: max(exp[h] / total, eps) for h in haps}
            s = sum(pi_new.values())
            pi_new = {h: pi_new[h] / s for h in haps}
        delta = max(abs(pi_new[h] - pi[h]) for h in haps)
        pi = pi_new
        if delta < tol:
            break
    return pi


def bootstrap_ci(rows: List[Tuple[str, Tuple[str, ...], int]], haps: List[str], boot: int, seed: int) -> Dict[str, Tuple[float, float]]:
    if boot <= 0:
        return {h: (math.nan, math.nan) for h in haps}
    rng = random.Random(seed)
    boots = []
    for _ in range(boot):
        b_rows = multinomial_resample(rows, rng)
        boots.append(em_pi_from_rows(b_rows, haps))
    ci = {}
    for h in haps:
        vals = sorted(b[h] for b in boots)
        lo = vals[int(0.025 * (len(vals) - 1))]
        hi = vals[int(0.975 * (len(vals) - 1))]
        ci[h] = (lo, hi)
    return ci


def write_pi(path: Path, pi: Dict[str, float], ci: Dict[str, Tuple[float, float]], total_reads: int, hap_labels: Dict[str, str]):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['hap_id', 'hap_label', 'pi', 'ci_low', 'ci_high', 'total_reads'])
        for h in pi.keys():
            lo, hi = ci.get(h, (math.nan, math.nan))
            w.writerow([h, hap_labels[h], f"{pi[h]:.6g}", f"{lo:.6g}", f"{hi:.6g}", total_reads])


def read_pi(path: Path) -> Tuple[Dict[str, Dict[str, float]], int]:
    data = {}
    total_reads = 0
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            h = row['hap_id']
            data[h] = {
                'pi': float(row['pi']),
                'ci_low': float(row['ci_low']),
                'ci_high': float(row['ci_high']),
            }
            total_reads = int(float(row['total_reads']))
    return data, total_reads


def summarize(cws_path: Path, pi_path: Path, panel: str, haps: List[str]) -> Dict[str, float]:
    rows = read_cws_rows(cws_path)
    total_reads = sum(c for _, _, c in rows)
    singleton = sum(c for _, hs, c in rows if len(hs) == 1)
    fullset = sum(c for _, hs, c in rows if len(hs) == len(haps))
    pi, _ = read_pi(pi_path)
    rec: Dict[str, float] = {
        'panel': panel,
        'total_reads': int(round(total_reads)),
        'singleton_frac': (singleton / total_reads) if total_reads else 0.0,
        'fullset_frac': (fullset / total_reads) if total_reads else 0.0,
    }
    ci_sum = 0.0
    for h in haps:
        rec[f'pi_{h}'] = pi[h]['pi']
        rec[f'w_{h}'] = pi[h]['ci_high'] - pi[h]['ci_low']
        ci_sum += rec[f'w_{h}']
    rec['ci_sum'] = ci_sum
    return rec


CI_TIE_MULT = 1.20


def panel_sort_key(rec: Dict[str, float]):
    # global ranking table for reporting: lower ci_sum first, then higher reads
    return (rec['ci_sum'], -rec['total_reads'], -rec['singleton_frac'], rec['fullset_frac'])


def tie_range_sort_key(rec: Dict[str, float]):
    # within acceptable ci range: prioritize higher reads first
    return (-rec['total_reads'], rec['ci_sum'], -rec['singleton_frac'], rec['fullset_frac'])


def open_maybe_gz(path: Path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')


def read_eq_rows(path: Path, min_count: float = 1.0) -> List[Tuple[str, Tuple[str, ...], int]]:
    with open_maybe_gz(path) as f:
        lines = [x.rstrip('\n') for x in f if x.strip() != '']
    if len(lines) < 3:
        raise ValueError(f'eq_classes too short: {path}')
    n_targets = int(lines[0])
    n_eq = int(lines[1])
    targets = lines[2:2+n_targets]
    eq_lines = lines[2+n_targets:2+n_targets+n_eq]
    parsed = []
    for t in targets:
        m = TARGET_RE.match(t)
        parsed.append((m.group(1), m.group(2)) if m else None)
    rows: List[Tuple[str, Tuple[str, ...], int]] = []
    for line in eq_lines:
        toks = line.split('\t')
        if len(toks) == 1:
            toks = line.split()
        if len(toks) < 2:
            continue
        k = int(toks[0])
        ids = list(map(int, toks[1:1+k]))
        cnt = float(toks[-1])
        if cnt < min_count:
            continue
        windows, haps = [], []
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
        rows.append(("+".join(sorted(set(windows))), tuple(sorted(set(haps))), int(round(cnt))))
    if not rows:
        raise ValueError(f'No usable eq rows after parsing/filtering: {path}')
    return rows


def aggregate_rows(rows: List[Tuple[str, Tuple[str, ...], int]]) -> List[Tuple[str, Tuple[str, ...], int]]:
    agg: Dict[Tuple[str, Tuple[str, ...]], int] = {}
    for w, hs, c in rows:
        if c > 0:
            agg[(w, hs)] = agg.get((w, hs), 0) + c
    return [(w, hs, c) for (w, hs), c in agg.items()]


def ec_bootstrap(eq_rows: List[Tuple[str, Tuple[str, ...], int]], haps: List[str], boot: int, seed: int):
    total_reads = sum(c for _, _, c in eq_rows)
    cws_rows = aggregate_rows(eq_rows)
    pi_hat = em_pi_from_rows(cws_rows, haps)
    ci = bootstrap_ci(eq_rows, haps, boot, seed)
    return cws_rows, pi_hat, ci, total_reads


def write_wide_tsv(path: Path, rows: List[Dict[str, object]], fields: List[str]):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main():
    args = parse_args()
    nt = resolve_nt(args.nt)
    print(f"Detected CPU cores: {os.cpu_count() or 1}")
    print(f"--nt {args.nt} resolved to: {nt}")
    outdir = Path(args.outdir).resolve()
    ensure_dirs(outdir)

    hap_rows = read_hap_map(outdir / 'build' / 'hap_map.tsv')
    haps = [hid for hid, _ in hap_rows]
    hap_labels = {hid: lab for hid, lab in hap_rows}

    panel_filter = set(parse_listish(args.panels) or []) or None
    sample_filter = set(parse_listish(args.samples) or []) or None
    panels = discover_panels(outdir, panel_filter)
    samples = discover_samples(outdir, panels[0], sample_filter)
    print(f"Panels: {','.join(panels)}")
    print(f"Samples: {','.join(samples)}")

    # Step 5: merge ALL3
    for p in panels:
        out = merge_all3(outdir / 'calib' / 'cws' / p, samples)
        print(f"[OK] Wrote {out}")

    # Step 6: EM with CwS bootstrap for screening
    tasks = []
    for p in panels:
        for s in samples:
            tasks.append((p, s, outdir / 'calib' / 'cws' / p / f'{s}.CwS.tsv', outdir / 'calib' / 'em' / f'{p}.{s}.pi.tsv', args.seed))
        tasks.append((p, 'ALL3', outdir / 'calib' / 'cws' / p / 'ALL3.CwS.tsv', outdir / 'calib' / 'em' / f'{p}.ALL3.pi.tsv', args.seed))

    def em_worker(task):
        panel, sample, cws_path, out_path, seed = task
        rows = read_cws_rows(cws_path)
        total_reads = sum(c for _, _, c in rows)
        pi = em_pi_from_rows(rows, haps)
        ci = bootstrap_ci(rows, haps, args.em_boot, seed)
        write_pi(out_path, pi, ci, total_reads, hap_labels)
        return panel, sample, out_path

    workers = max(1, min(nt, len(tasks)))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(em_worker, t) for t in tasks]
        for fu in as_completed(futs):
            panel, sample, out_path = fu.result()
            print(f"[OK] EM {panel} {sample}: {out_path}")

    # Step 7: panel summaries and selection
    panel_rows = []
    for p in panels:
        panel_rows.append(summarize(outdir / 'calib' / 'cws' / p / 'ALL3.CwS.tsv', outdir / 'calib' / 'em' / f'{p}.ALL3.pi.tsv', p, haps))
    ranked = sorted(panel_rows, key=panel_sort_key)

    ci_best = min((rec['ci_sum'] for rec in ranked), default=math.nan)
    ci_cutoff = ci_best * CI_TIE_MULT if ranked else math.nan
    near_tie = [rec for rec in ranked if rec['ci_sum'] <= ci_cutoff] if ranked else []
    near_tie_sorted = sorted(near_tie, key=tie_range_sort_key)

    for i, rec in enumerate(ranked, start=1):
        rec['rank'] = i
        rec['selected'] = 'no'
        rec['recommendation'] = ''
        rec['status'] = 'OK'

    best_panel = ''
    runner_up = ''
    selected_panels = []

    if args.panel:
        best_panel = args.panel
        selected_panels = [args.panel]
        for rec in ranked:
            if rec['panel'] == args.panel:
                rec['selected'] = 'yes'
                rec['recommendation'] = 'manual_override'
                break
        remaining = [rec for rec in ranked if rec['panel'] != args.panel]
        if remaining:
            runner_up = remaining[0]['panel']
            remaining[0]['recommendation'] = 'runner_up'
    else:
        if near_tie_sorted:
            best_panel = near_tie_sorted[0]['panel']
            selected_panels = [rec['panel'] for rec in near_tie_sorted[:max(1, args.ecboot_top)]]
            if len(near_tie_sorted) > 1:
                runner_up = near_tie_sorted[1]['panel']
            elif len(ranked) > 1:
                runner_up = ranked[1]['panel']

            for rec in ranked:
                if rec['panel'] == best_panel:
                    rec['selected'] = 'yes'
                    rec['recommendation'] = 'best_panel'
                elif rec['panel'] == runner_up:
                    rec['recommendation'] = 'runner_up'
                elif rec['ci_sum'] <= ci_cutoff:
                    rec['recommendation'] = 'within_1.20x_ci_best'
        elif ranked:
            best_panel = ranked[0]['panel']
            selected_panels = [best_panel]
            ranked[0]['selected'] = 'yes'
            ranked[0]['recommendation'] = 'best_panel'
            if len(ranked) > 1:
                runner_up = ranked[1]['panel']
                ranked[1]['recommendation'] = 'runner_up'

    panel_fields = ['panel', 'total_reads', 'singleton_frac', 'fullset_frac', 'ci_sum', 'rank', 'selected', 'recommendation', 'status']
    for h in haps:
        panel_fields.extend([f'pi_{h}', f'w_{h}'])
    panel_summary_path = outdir / 'calib' / 'summary' / 'panel_summary.tsv'
    write_wide_tsv(panel_summary_path, ranked, panel_fields)

    sample_rows = []
    for p in panels:
        for s in samples:
            rec = summarize(outdir / 'calib' / 'cws' / p / f'{s}.CwS.tsv', outdir / 'calib' / 'em' / f'{p}.{s}.pi.tsv', p, haps)
            rec['sample'] = s
            sample_rows.append(rec)
    sample_fields = ['panel', 'sample', 'total_reads', 'singleton_frac', 'fullset_frac', 'ci_sum']
    for h in haps:
        sample_fields.extend([f'pi_{h}', f'w_{h}'])
    sample_summary_path = outdir / 'calib' / 'summary' / 'sample_summary.tsv'
    write_wide_tsv(sample_summary_path, sorted(sample_rows, key=lambda x: (x['panel'], x['sample'])), sample_fields)


    final_summary_path = outdir / 'calib' / 'summary' / 'final_summary.tsv'
    with open(final_summary_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['best_panel', 'runner_up', 'stop_after', 'selected_panels'])
        w.writerow([best_panel, runner_up, args.stop_after, ','.join(selected_panels)])

    if args.stop_after == 'select':
        print("Stopped after panel selection.")
        print(f"Recommended best panel: {best_panel}")
        print(f"Runner-up: {runner_up}")
        print(f"Panel summary: {panel_summary_path}")
        print(f"Sample summary: {sample_summary_path}")
        print(f"Final summary: {final_summary_path}")
        return

    # Step 8: EC bootstrap on selected panel(s)
    ec_tasks = []
    for p in selected_panels:
        if not p:
            continue
        for s in samples:
            eqs = sorted((outdir / 'calib' / 'quant' / p / s / 'aux_info').glob('eq_classes.txt*'))
            if not eqs:
                raise FileNotFoundError(f"Missing eq_classes for {p} {s}")
            ec_tasks.append((p, s, eqs, outdir / 'calib' / 'ecboot' / f'{p}.{s}', args.seed))
        all_eqs = []
        for s in samples:
            eqs = sorted((outdir / 'calib' / 'quant' / p / s / 'aux_info').glob('eq_classes.txt*'))
            if not eqs:
                raise FileNotFoundError(f"Missing eq_classes for {p} {s}")
            all_eqs.extend(eqs)
        ec_tasks.append((p, 'ALL3', all_eqs, outdir / 'calib' / 'ecboot' / f'{p}.ALL3', args.seed))

    def ec_worker(task):
        panel, sample, eq_paths, out_prefix, seed = task
        eq_rows = []
        for ep in eq_paths:
            eq_rows.extend(read_eq_rows(ep))
        cws_rows, pi, ci, total_reads = ec_bootstrap(eq_rows, haps, args.ec_boot, seed)
        write_cws_rows(Path(str(out_prefix) + '.cws.tsv'), cws_rows)
        write_pi(Path(str(out_prefix) + '.pi.tsv'), pi, ci, total_reads, hap_labels)
        return panel, sample, str(out_prefix) + '.pi.tsv'

    workers = max(1, min(nt, len(ec_tasks)))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(ec_worker, t) for t in ec_tasks]
        for fu in as_completed(futs):
            panel, sample, out_path = fu.result()
            print(f"[OK] EC bootstrap {panel} {sample}: {out_path}")

    final_haps_path = outdir / 'calib' / 'summary' / 'final_hap_fractions.tsv'
    final_rows = []
    for p in selected_panels:
        for tag in samples + ['ALL3']:
            pi_path = outdir / 'calib' / 'ecboot' / (f'{p}.{tag}.pi.tsv' if tag != 'ALL3' else f'{p}.ALL3.pi.tsv')
            pi, total_reads = read_pi(pi_path)
            for h in haps:
                final_rows.append({
                    'panel': p,
                    'sample': tag,
                    'hap_id': h,
                    'hap_label': hap_labels[h],
                    'pi': pi[h]['pi'],
                    'ci_low': pi[h]['ci_low'],
                    'ci_high': pi[h]['ci_high'],
                    'total_reads': total_reads,
                })
    with open(final_haps_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['panel','sample','hap_id','hap_label','pi','ci_low','ci_high','total_reads'], delimiter='\t')
        w.writeheader()
        for row in final_rows:
            w.writerow(row)

    print("Inference pipeline finished.")
    print(f"Panel summary: {panel_summary_path}")
    print(f"Sample summary: {sample_summary_path}")
    print(f"Final summary: {final_summary_path}")
    print(f"Final hap fractions: {final_haps_path}")
    print(f"Best panel: {best_panel}")
    print(f"EC bootstrap directory: {outdir / 'calib' / 'ecboot'}")


if __name__ == '__main__':
    main()
