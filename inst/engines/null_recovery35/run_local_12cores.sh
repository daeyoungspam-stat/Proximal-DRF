#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PDR_NULL_ROOT="$ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PDR_PROGRESS=1
if [[ -z "${R_LIBS_USER:-}" ]]; then
  LOCAL_LIB="$ROOT/../../Derivative R study/joint_simulation/.Rlib"
  if [[ -d "$LOCAL_LIB" ]]; then export R_LIBS_USER="$LOCAL_LIB"; fi
fi
cd "$ROOT"
RSCRIPT="${RSCRIPT:-Rscript}"
if [[ "${1:-}" == "--detach" ]]; then
  shift
  mkdir -p "$ROOT/launch_logs"
  LOG="$ROOT/launch_logs/local_$(date -u +%Y%m%dT%H%M%SZ).log"
  # A new process session also survives callers that clean up their process
  # group when the launching shell exits. nohup alone does not provide this.
  python3 - "$ROOT" "$LOG" "$@" <<'PY'
import shutil
import subprocess
import sys

root, log_path, *arguments = sys.argv[1:]
command = ["bash", root + "/run_local_12cores.sh", *arguments]
if shutil.which("caffeinate"):
    command = ["caffeinate", "-i", *command]
with open(log_path, "ab", buffering=0) as log:
    process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log,
                               stderr=subprocess.STDOUT, start_new_session=True)
print(f"Launched PID {process.pid}. Log: {log_path}")
PY
  exit 0
fi
"$RSCRIPT" "$ROOT/scripts/run_local.R" --workers 12 "$@"
"$RSCRIPT" "$ROOT/scripts/combine_results.R" "$@"
"$RSCRIPT" "$ROOT/scripts/make_report.R" "$@"
