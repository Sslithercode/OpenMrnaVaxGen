# ── Melanoma mRNA Vaccine Pipeline ────────────────────────────────────────────
# Single image containing all dependencies for Steps 1–7 (and Step 8).
#
# Base: CUDA 12.4.1 + cuDNN 9 on Ubuntu 22.04
#   • Steps 4–6: MHCflurry (TensorFlow) — GPU accelerated
#   • Step 7:    VaxPress — GPU accelerated via ViennaRNA/LinearFold
#
# Build:
#   docker build -t melanoma-pipeline .
#
# Launch Streamlit UI (default):
#   docker compose up app
#   open http://localhost:8501
#
# Run a step via CLI:
#   docker compose run --rm pipeline python3 scripts/preprocess.py
#   docker compose run --rm pipeline python3 scripts/variant.py
#
# Step 3 OptiType is a separate service — see docker-compose.yml.
# ──────────────────────────────────────────────────────────────────────────────

FROM nvidia/cuda:12.4.1-cudnn9-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# ── System packages ────────────────────────────────────────────────────────────
# deadsnakes PPA: Python 3.12
# universe repo:  vienna-rna, samtools
# manual:         Java 17 (for GATK)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        curl \
        wget \
        gnupg \
        ca-certificates \
        unzip \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && add-apt-repository -y universe \
    && apt-get update && \
    apt-get install -y --no-install-recommends \
        python3.12 \
        python3.12-dev \
        python3.12-distutils \
        openjdk-17-jre-headless \
        samtools \
        tabix \
        bcftools \
        vienna-rna \
        docker.io \
        bwa \
        fastqc \
        sra-toolkit \
        bedtools \
        kallisto \
        pigz \
    && rm -rf /var/lib/apt/lists/*

# ── Python setup ───────────────────────────────────────────────────────────────
RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.12 && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python3 && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python

# ── GATK 4.5.0.0 ──────────────────────────────────────────────────────────────
# Downloaded and placed at the path the scripts expect via Path.home() resolution
# (scripts use Path.home() / "melanoma-pipeline/tools/gatk-4.5.0.0/...")
RUN mkdir -p /root/melanoma-pipeline/tools && \
    wget -q \
        https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip \
        -O /tmp/gatk.zip && \
    unzip -q /tmp/gatk.zip -d /root/melanoma-pipeline/tools && \
    rm /tmp/gatk.zip

# ── Python dependencies ────────────────────────────────────────────────────────
# Install in layers ordered by size/stability so Docker cache is maximised:
#   1. TensorFlow (largest, least likely to change version-to-version)
#   2. PyTorch (for Step 8 CodonFM)
#   3. MHCflurry + bioinformatics libs
#   4. VaxPress + mRNA design libs
#   5. Pipeline utilities

RUN python3 -m pip install --no-cache-dir \
    tensorflow==2.20.0

RUN python3 -m pip install --no-cache-dir \
    torch --index-url https://download.pytorch.org/whl/cu124

RUN python3 -m pip install --no-cache-dir \
    mhcflurry>=2.1.5 \
    pysam>=0.22.0 \
    biopython>=1.81

RUN python3 -m pip install --no-cache-dir \
    vaxpress>=0.9 \
    transformers>=5.3.0 \
    sentencepiece>=0.2.1

RUN python3 -m pip install --no-cache-dir \
    pandas>=2.0.0 \
    numpy>=1.24.0 \
    matplotlib>=3.7.0 \
    seaborn>=0.12.0 \
    requests>=2.31.0 \
    tqdm>=4.65.0 \
    reportlab>=4.0.0 \
    streamlit>=1.35.0 \
    linearfold-unofficial>=0.1

# ── Download MHCflurry pretrained models ───────────────────────────────────────
# Baked into the image so containers are fully self-contained at runtime.
RUN mhcflurry-downloads fetch

# ── Copy pipeline source ───────────────────────────────────────────────────────
WORKDIR /root/melanoma-pipeline
COPY scripts/    scripts/
COPY src/        src/
COPY app.py      app.py
COPY pyproject.toml pyproject.toml

# ── Streamlit port ─────────────────────────────────────────────────────────────
EXPOSE 8501

# ── Default: launch Streamlit UI ───────────────────────────────────────────────
CMD ["streamlit", "run", "app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true"]
