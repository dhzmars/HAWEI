#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT=""
OUTDIR=""
FQ_DIR=""
FQ_MANIFEST=""
LENS="300,400,800,1200"
NT="auto"
STOP_AFTER="ecboot"
PANELS=""
SAMPLES=""
MIN_UNGAPPED="80"
SALMON_BIN="salmon"
KMER="31"

usage() {
  cat <<EOF
Usage:
  bash HAWEI.sh -i <align.fasta> -o <route_dir> (--fq <fastq_dir> | --fq-manifest <samples.tsv>) [options]

Required:
  -i, --input           aligned full-length hap FASTA
  -o, --outdir          output route directory
  --fq                  directory containing sample_1.fastq.gz and sample_2.fastq.gz
  --fq-manifest         TSV manifest with columns: sample r1 r2

Optional:
  --len                 window lengths [default: 300,400,800,1200]
  --nt                  threads, supports auto/8/12 [default: auto]
  --stop-after          build|cws|select|ecboot [default: ecboot]
  --panels              optional panel list for downstream steps
  --samples             optional sample list
  --min-ungapped        minimum ungapped nt per hap-window [default: 80]
  --salmon              salmon executable [default: salmon]
  --kmer                salmon index k-mer [default: 31]
  -h, --help            show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INPUT="$2"; shift 2 ;;
    -o|--outdir)
      OUTDIR="$2"; shift 2 ;;
    --fq)
      FQ_DIR="$2"; shift 2 ;;
    --fq-manifest)
      FQ_MANIFEST="$2"; shift 2 ;;
    --len)
      LENS="$2"; shift 2 ;;
    --nt)
      NT="$2"; shift 2 ;;
    --stop-after)
      STOP_AFTER="$2"; shift 2 ;;
    --panels)
      PANELS="$2"; shift 2 ;;
    --samples)
      SAMPLES="$2"; shift 2 ;;
    --min-ungapped)
      MIN_UNGAPPED="$2"; shift 2 ;;
    --salmon)
      SALMON_BIN="$2"; shift 2 ;;
    --kmer)
      KMER="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1 ;;
  esac
done

if [[ -z "$INPUT" || -z "$OUTDIR" ]]; then
  echo "Error: -i/--input and -o/--outdir are required." >&2
  usage >&2
  exit 1
fi

if [[ -n "$FQ_DIR" && -n "$FQ_MANIFEST" ]]; then
  echo "Error: use only one of --fq or --fq-manifest." >&2
  exit 1
fi

if [[ -z "$FQ_DIR" && -z "$FQ_MANIFEST" ]]; then
  echo "Error: one of --fq or --fq-manifest is required." >&2
  exit 1
fi

case "$STOP_AFTER" in
  build|cws|select|ecboot) ;;
  *) echo "Error: --stop-after must be one of build|cws|select|ecboot" >&2; exit 1 ;;
esac

INPUT="$(readlink -f "$INPUT")"
OUTDIR="$(readlink -m "$OUTDIR")"
if [[ -n "$FQ_DIR" ]]; then FQ_DIR="$(readlink -f "$FQ_DIR")"; fi
if [[ -n "$FQ_MANIFEST" ]]; then FQ_MANIFEST="$(readlink -f "$FQ_MANIFEST")"; fi

mkdir -p "$OUTDIR/logs"
TS="$(date +%Y%m%d_%H%M%S)"
MAIN_LOG="$OUTDIR/logs/HAWEI.${TS}.log"
BUILD_LOG="$OUTDIR/logs/build.log"
MAINPANELS_LOG="$OUTDIR/logs/main_panels.log"
INFER_LOG="$OUTDIR/logs/inference.log"

log_msg() {
  echo "$1" | tee -a "$MAIN_LOG"
}

run_stage() {
  local stage="$1"
  local stage_log="$2"
  shift 2
  log_msg "[STAGE] ${stage}"
  log_msg "[CMD] $*"
  set +e
  "$@" 2>&1 | tee "$stage_log" | tee -a "$MAIN_LOG"
  local rc=${PIPESTATUS[0]}
  set -e
  if [[ $rc -ne 0 ]]; then
    log_msg "[FAIL] ${stage} exit code ${rc}"
    exit $rc
  fi
  log_msg "[OK] ${stage}"
}

log_msg "HAWEI started"
log_msg "Input alignment: $INPUT"
log_msg "Outdir: $OUTDIR"
if [[ -n "$FQ_DIR" ]]; then log_msg "FASTQ directory: $FQ_DIR"; fi
if [[ -n "$FQ_MANIFEST" ]]; then log_msg "FASTQ manifest: $FQ_MANIFEST"; fi
log_msg "Window lengths: $LENS"
log_msg "--nt: $NT"
log_msg "--stop-after: $STOP_AFTER"

# Stage 1: build
run_stage "build" "$BUILD_LOG" \
  python3 "$SCRIPT_DIR/make_multiL_windows_from_align.py" \
  -i "$INPUT" -o "$OUTDIR" --len "$LENS" --min-ungapped "$MIN_UNGAPPED" --nt "$NT"

if [[ "$STOP_AFTER" == "build" ]]; then
  log_msg "Stopped after build"
  echo "Stopped after build."
  echo "Main log: $MAIN_LOG"
  echo "Build windows: $OUTDIR/build/windows_multiL.clean.fasta"
  exit 0
fi

# Stage 2: main panels
MAIN_ARGS=(python3 "$SCRIPT_DIR/run_main_panels.py" -o "$OUTDIR" --nt "$NT" --salmon "$SALMON_BIN" --kmer "$KMER")
if [[ -n "$FQ_DIR" ]]; then
  MAIN_ARGS+=(--fq "$FQ_DIR")
else
  MAIN_ARGS+=(--fq-manifest "$FQ_MANIFEST")
fi
if [[ -n "$PANELS" ]]; then MAIN_ARGS+=(--panels "$PANELS"); fi
if [[ -n "$SAMPLES" ]]; then MAIN_ARGS+=(--samples "$SAMPLES"); fi
run_stage "main_panels" "$MAINPANELS_LOG" "${MAIN_ARGS[@]}"

if [[ "$STOP_AFTER" == "cws" ]]; then
  log_msg "Stopped after cws"
  echo "Stopped after cws."
  echo "Main log: $MAIN_LOG"
  echo "CwS outputs: $OUTDIR/calib/cws/"
  echo "Run summary: $OUTDIR/calib/run_main_panels.summary.tsv"
  exit 0
fi

# Stage 3: inference
INFER_ARGS=(python3 "$SCRIPT_DIR/run_panel_inference.py" -o "$OUTDIR" --nt "$NT")
if [[ -n "$PANELS" ]]; then INFER_ARGS+=(--panels "$PANELS"); fi
if [[ -n "$SAMPLES" ]]; then INFER_ARGS+=(--samples "$SAMPLES"); fi
if [[ "$STOP_AFTER" == "select" ]]; then
  INFER_ARGS+=(--stop-after select)
fi
run_stage "inference" "$INFER_LOG" "${INFER_ARGS[@]}"

BEST_PANEL=""
FINAL_SUMMARY="$OUTDIR/calib/summary/final_summary.tsv"
if [[ -f "$FINAL_SUMMARY" ]]; then
  BEST_PANEL="$(tail -n +2 "$FINAL_SUMMARY" | head -n 1 | cut -f1 || true)"
fi

log_msg "Pipeline finished"

echo "Pipeline finished."
echo "Outdir: $OUTDIR"
echo "Main log: $MAIN_LOG"
echo "Build windows: $OUTDIR/build/windows_multiL.clean.fasta"
echo "Hap map: $OUTDIR/build/hap_map.tsv"
echo "Main panel FASTA: $OUTDIR/calib/fasta/"
echo "Quant outputs: $OUTDIR/calib/quant/"
echo "CwS outputs: $OUTDIR/calib/cws/"
echo "EM outputs: $OUTDIR/calib/em/"
echo "EC bootstrap outputs: $OUTDIR/calib/ecboot/"
echo "Panel summary: $OUTDIR/calib/summary/panel_summary.tsv"
echo "Sample summary: $OUTDIR/calib/summary/sample_summary.tsv"
echo "Final summary: $OUTDIR/calib/summary/final_summary.tsv"
echo "Final hap fractions: $OUTDIR/calib/summary/final_hap_fractions.tsv"
if [[ -n "$BEST_PANEL" ]]; then
  echo "Best panel: $BEST_PANEL"
fi
