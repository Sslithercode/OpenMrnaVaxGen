"""
Download HCC1143 test data from HuggingFace.

Reads HF_TOKEN from .env in the project root.
Run from the project root:
    uv run --with huggingface_hub --with python-dotenv --no-project python scripts/download_data.py
"""

import os
from pathlib import Path

from dotenv import load_dotenv
from huggingface_hub import snapshot_download

ROOT = Path(__file__).resolve().parent.parent
load_dotenv(ROOT / ".env")

token = os.getenv("HF_TOKEN")
if not token:
    print("Warning: HF_TOKEN not set in .env — proceeding unauthenticated (rate limits apply)")

snapshot_download(
    repo_id="SlitherCode/hcc1143_cancer_and_normal_data",
    repo_type="dataset",
    allow_patterns="test/*",
    local_dir=str(ROOT / "data"),
    token=token or None,
)

print("Done — files saved to data/test/")
