"""
Download and decompress the hg38 reference genome from UCSC.

Output: reference/hg38.fa (~3.5 GB uncompressed)

Indexing (run after Docker image is built):
    docker compose run --rm index-reference

Run from the project root:
    uv run --with requests --with tqdm --no-project python scripts/download_reference.py
"""

import gzip
import shutil
from pathlib import Path

import requests
from tqdm import tqdm

URL = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
ROOT = Path(__file__).resolve().parent.parent
OUT_GZ = ROOT / "reference" / "hg38.fa.gz"
OUT_FA = ROOT / "reference" / "hg38.fa"


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(dest, "wb") as f, tqdm(
            desc="Downloading hg38.fa.gz",
            total=total,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)
                bar.update(len(chunk))


def decompress(src: Path, dest: Path) -> None:
    print(f"Decompressing {src.name} ...")
    with gzip.open(src, "rb") as f_in, open(dest, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    src.unlink()
    print(f"Done — {dest}")


if __name__ == "__main__":
    if OUT_FA.exists():
        print(f"{OUT_FA} already exists, skipping.")
    else:
        if not OUT_GZ.exists():
            download(URL, OUT_GZ)
        decompress(OUT_GZ, OUT_FA)
        print("\nNext step — index inside Docker:")
        print("  docker compose run --rm index-reference")
