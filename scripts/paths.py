"""
Central path resolver for the melanoma-pipeline.

Every step script imports from here instead of hardcoding its own
results/stepN paths.  The run directory is resolved in this priority order:

  1. PIPELINE_RUN env var          — explicit override / parallel runs
  2. results/.current_run file     — same-session continuity (no env var needed)
  3. Auto-generate a new run ID    — first step of a fresh run

Format:  run_{SAMPLE_ID}_{YYYYMMDD}_{HHMMSS}
Example: run_HCC1143_20260318_143022

Usage in step scripts:
    from paths import step_dir, RUN_DIR, SAMPLE_ID
    OUT_DIR      = step_dir(5)
    STEP4_OUTPUT = step_dir(4) / "candidate_neoantigens.tsv"

To force a brand-new run (call only from step 1 / main.py):
    from paths import new_run
    new_run()
"""

import os
from datetime import datetime
from pathlib import Path

# ── Base locations ─────────────────────────────────────────────────────────────

BASE    = Path.home() / "melanoma-pipeline"
RESULTS = BASE / "results"

# Sample / run identity — override via env vars for non-HCC1143 samples
SAMPLE_ID    = os.environ.get("PIPELINE_SAMPLE",  "HCC1143")
TUMOR_TYPE   = os.environ.get("PIPELINE_TUMOR",   "Triple-negative breast cancer (TNBC)")
PIPELINE_VER = "1.0.0"

# Sentinel file that persists the active run ID across individual step invocations
_RUN_FILE = RESULTS / ".current_run"


# ── Run ID resolution ──────────────────────────────────────────────────────────

def _generate_run_id() -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"run_{SAMPLE_ID}_{ts}"


def new_run() -> str:
    """
    Generate a fresh run ID and write it to the sentinel file.
    Call this only from step 1 or main.py — not from downstream steps.
    Returns the new run ID.
    """
    run_id = _generate_run_id()
    RESULTS.mkdir(parents=True, exist_ok=True)
    _RUN_FILE.write_text(run_id + "\n")
    print(f"[paths] New run: {run_id}")
    return run_id


def _resolve_run_id() -> str:
    # 1. Explicit env var (parallel runs, CI, reproducibility)
    env = os.environ.get("PIPELINE_RUN", "").strip()
    if env:
        return env

    # 2. Sentinel file written by the most recent call to new_run()
    if _RUN_FILE.exists():
        run_id = _RUN_FILE.read_text().strip()
        if run_id:
            return run_id

    # 3. No active run — generate one (allows running steps standalone)
    return new_run()


# ── Public API ─────────────────────────────────────────────────────────────────

RUN_ID  = _resolve_run_id()
RUN_DIR = RESULTS / RUN_ID


def step_dir(n: int) -> Path:
    """Return the output directory for pipeline step n within the active run."""
    return RUN_DIR / f"step{n}"
