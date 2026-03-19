# syntax=docker/dockerfile:1
# ── STAGE 1: Builder ──────────────────────────────────────────────────────────
FROM python:3.12-slim-bookworm AS builder

# Install build-only dependencies
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
    wget unzip git build-essential curl

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /app

# 1. Download GATK 4.5
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip -O /tmp/gatk.zip \
    && unzip -q /tmp/gatk.zip -d /app/tools && rm /tmp/gatk.zip

# 2. Clone LinearFold
RUN git clone --depth 1 https://github.com/LinearFold/LinearFold.git /opt/LinearFold

# ── STAGE 2: Final Runtime ────────────────────────────────────────────────────
FROM python:3.12-slim-bookworm

# Environment setup
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/app/.venv/bin:$PATH" \
    DEBIAN_FRONTEND=noninteractive \
    UV_COMPILE_BYTECODE=1 \
    UV_PYTHON_DOWNLOADS=never \
    UV_PROJECT_ENVIRONMENT=/app/.venv

WORKDIR /app

# Copy uv for runtime use (entrypoint runs uv sync against mounted venv volume)
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Install Runtime Bioinformatics Tools
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless \
    libgomp1 \
    libjemalloc2 \
    samtools \
    tabix \
    bcftools \
    bwa \
    bedtools \
    pigz \
    curl

# Copy pre-built assets from builder
COPY --from=builder /app/tools /app/tools
COPY --from=builder /opt/LinearFold /LinearFold

# Copy lockfile so uv sync can run at container startup
COPY pyproject.toml uv.lock ./

# Copy Application Code
COPY scripts/ scripts/
COPY src/ src/
COPY app.py .

# Entrypoint: syncs venv + downloads models on first run, fast on subsequent runs
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

EXPOSE 8501

ENTRYPOINT ["/entrypoint.sh"]
CMD ["streamlit", "run", "app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true"]
