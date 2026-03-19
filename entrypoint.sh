#!/bin/bash
set -e

# Sync packages into the persistent venv volume.
# Fast if uv.lock is unchanged; installs only the diff otherwise.
echo "--- Syncing packages ---"
uv sync --frozen --no-install-project --no-dev

# MHCflurry models — download if missing (idempotent, skips if already present)
if [ ! -d "/root/.cache/mhcflurry/4" ] || [ -z "$(ls -A /root/.cache/mhcflurry/4 2>/dev/null)" ]; then
    echo "--- Downloading MHCflurry models ---"
    /app/.venv/bin/mhcflurry-downloads fetch
fi

# pyensembl GRCh38 release 109 — defaults to /root/.cache/pyensembl
if [ ! -d "/root/.cache/pyensembl/GRCh38/ensembl109" ]; then
    echo "--- Installing pyensembl GRCh38 r109 (first run) ---"
    /app/.venv/bin/pyensembl install --release 109 --species homo_sapiens
fi

exec "$@"
