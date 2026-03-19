# ── STAGE 1: Builder ──────────────────────────────────────────────────────────
FROM python:3.12-slim-bookworm AS builder

# Install build-only dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget unzip git build-essential curl \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /app

# 1. Download GATK 4.5
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip -O /tmp/gatk.zip \
    && unzip -q /tmp/gatk.zip -d /app/tools && rm /tmp/gatk.zip

# 2. Clone LinearFold
RUN git clone --depth 1 https://github.com/LinearFold/LinearFold.git /opt/LinearFold

# 3. Install Python deps into a virtualenv at /app/.venv so shebangs match runtime path
ENV UV_COMPILE_BYTECODE=1 \
    UV_PYTHON_DOWNLOADS=never
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --no-install-project --no-dev

# ── STAGE 2: Final Runtime ────────────────────────────────────────────────────
FROM python:3.12-slim-bookworm

# Environment setup
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/app/.venv/bin:$PATH" \
    DEBIAN_FRONTEND=noninteractive

WORKDIR /app

# Copy uv for runtime utility
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Install Runtime Bioinformatics Tools
# Note: Debian 'bookworm' (the base of this slim image) includes most 
# bio-tools in its default main/universe-equivalent repos.
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless \
    samtools \
    tabix \
    bcftools \
    bwa \
    bedtools \
    pigz \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy pre-built assets from builder
COPY --from=builder /app/tools /app/tools
COPY --from=builder /opt/LinearFold /LinearFold
COPY --from=builder /app/.venv /app/.venv

# Copy Application Code
COPY scripts/ scripts/
COPY src/ src/
COPY app.py .

# Download MHCflurry models (layer this last as it's large)
RUN /app/.venv/bin/mhcflurry-downloads fetch

EXPOSE 8501

CMD ["streamlit", "run", "app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true"]