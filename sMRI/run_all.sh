#!/bin/bash
set -euo pipefail

# ============================================================
# Master pipeline runner
# Runs S1 -> S5 in order and prints progress messages.
# Each step writes stdout/stderr to a timestamped log file.
# If any step fails, the script stops and prints the last lines
# of the corresponding log for quick debugging.
# ============================================================

# NOTE:
# Some conda activate/deactivate hooks may reference variables that
# are unset. If your sub-scripts use "set -u" and conda, you can
# temporarily disable nounset around conda calls inside those scripts.
# This wrapper uses bash -lc to run each command in a login shell.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Directory to store logs from all steps
LOG_DIR="$SCRIPT_DIR/logs_run_all"
mkdir -p "$LOG_DIR"

# Timestamp used in log filenames
ts="$(date +%Y%m%d_%H%M%S)"

run_step () {
  local step_id="$1"     # Step tag, e.g., S1, S2...
  local name="$2"        # Human-readable step name
  local cmd="$3"         # Command to execute
  local log="$LOG_DIR/${ts}_${step_id}_${name}.log"

  echo "============================================================"
  echo "[$(date '+%F %T')] START  Step ${step_id}: ${name}"
  echo "Command: $cmd"
  echo "Log: $log"
  echo "------------------------------------------------------------"

  # Run the command and capture both stdout and stderr to the log file
  if bash -lc "$cmd" >"$log" 2>&1; then
    echo "[$(date '+%F %T')] OK     Step ${step_id}: ${name}"
  else
    echo "[$(date '+%F %T')] FAILED Step ${step_id}: ${name}"
    echo "---- Last 60 lines of log ----"
    tail -n 60 "$log" || true
    echo "--------------------------------"
    echo "Stopped at Step ${step_id}. Please check: $log"
    exit 1
  fi
}

# ============================================================
# Execute steps in order
# ============================================================
run_step "S1" "create_AAL_annot"               "bash S1_create_AAL_annot.sh"
run_step "S2" "ROI_analysis"                   "bash S2_ROI_analysis.sh"
run_step "S3" "AAL_olfactory_extract_simplify" "bash S3_AAL_olfactory_extract_simplify.sh"
run_step "S4" "extract_freesurfer_value"       "bash S4_extract_freesurfer_value.sh"
run_step "S5" "statistics_arrange"             "python S5_statistics_arrange.py"

echo "============================================================"
echo "ALL DONE ✅  Logs are in: $LOG_DIR"
echo "============================================================"
