# mRNA Vaccine Pipeline

> **Proof of concept / work in progress.** Not for clinical use.

An open-source computational pipeline for designing personalized mRNA vaccines for melanoma. Given tumor/normal sequencing data and a patient's HLA profile, it identifies tumor-specific neoantigens, ranks them by immunogenicity, and encodes the top candidates into an optimized mRNA vaccine sequence.

The UI is a [Streamlit](https://streamlit.io) app that lets you run and configure each step interactively. Docker is the recommended way to run it.

---

## Pipeline

```
Input: Tumor WES + Normal WES
       ↓
Step 1: Preprocessing          → Clean BAMs             [GATK MarkDuplicates + SortSam]
Step 2: Variant Calling        → Somatic VCF            [GATK Mutect2 + FilterMutectCalls]
Step 3: HLA Typing             → HLA alleles            [OptiType]
Step 4: Neoantigen Prediction  → Candidate peptides     [MHCflurry]
Step 5: Immunogenicity Ranking → Ranked neoantigen list [composite score]
Step 6: Epitope Ordering       → Ordered epitope string [TSP / greedy]
Step 7: mRNA Design            → Full mRNA sequence     [VaxPress + RNAfold]
Step 9: Report                 → vaccine_report.md
       ↓
Output: results/run_<sample>_<timestamp>/
```

**Entry A — full:** Raw FASTQs → Steps 1 → 2 → 3 → 4 → 5 → 6 → 7 → 9

**Entry B — pre-called:** MAF/VCF → Steps 3 → 4 → 5 → 6 → 7 → 9

---

## Prerequisites

The following are gitignored and must be set up manually before running.

### 1. Create required directories

```bash
mkdir -p tools data/test reference
```

### 2. Download GATK 4.5.0.0

```bash
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip -d tools/
rm gatk-4.5.0.0.zip
```

### 3. Add reference genomes

Drop the following into `reference/`:

| File | Used by |
|------|---------|
| `hg38.fa` + `hg38.fa.fai` + `hg38.dict` | Steps 2, 4 (GATK Mutect2 / MHCflurry) |
| `b37.20.21.fasta` + index | Step 1 (preprocessing) |

The hg38 reference is available from GATK's Google Cloud bucket (`gs://genomics-public-data/resources/broad/hg38/`) or from UCSC. Both require indexing with `samtools faidx` and `gatk CreateSequenceDictionary`.

### 4. Add input data

Drop your tumor and normal BAM/BAI files into `data/test/`:

```
data/test/
├── test_tumor.bam
├── test_tumor.bam.bai
├── normal_chr17.bam
└── normal_chr17.bam.bai
```

The pipeline defaults to these paths in the UI. You can override them per-step.

---

## Running with Docker (recommended)

### Launch the Streamlit UI

```bash
docker compose up app
```

Then open [http://localhost:8501](http://localhost:8501).

Volumes are mounted automatically:

| Host path | Container path |
|-----------|----------------|
| `./data` | `/root/melanoma-pipeline/data` |
| `./reference` | `/root/melanoma-pipeline/reference` |
| `./results` | `/root/melanoma-pipeline/results` |
| `./tools` | `/root/melanoma-pipeline/tools` |

Override paths with env vars:

```bash
DATA_DIR=/path/to/data REFERENCE_DIR=/path/to/ref docker compose up app
```

### Run individual steps via CLI

```bash
docker compose run --rm pipeline python3 scripts/preprocess.py
docker compose run --rm pipeline python3 scripts/variant.py
docker compose run --rm pipeline python3 scripts/neoantigen_prediction.py
docker compose run --rm pipeline python3 scripts/candidate_ranking.py
docker compose run --rm pipeline python3 scripts/epitope_ordering.py
docker compose run --rm pipeline python3 scripts/mrna_design.py
```

### Step 3: HLA Typing (three phases)

HLA typing requires OptiType, which runs as a separate container between two pipeline phases.

**Phase 1** — extract HLA reads from the normal BAM:
```bash
docker compose run --rm pipeline python3 scripts/hla_typing.py --extract
```

**Phase 2** — run OptiType (set `OPTITYPE_DATA_DIR` to the active run's `step3` dir):
```bash
OPTITYPE_DATA_DIR=./results/run_<id>/step3 \
OPTITYPE_SAMPLE=hcc1143_normal \
docker compose run --rm optitype
```

**Phase 3** — parse OptiType results:
```bash
docker compose run --rm pipeline python3 scripts/hla_typing.py --parse
```

Alternatively, if you already have HLA alleles, use the **Manual HLA entry** panel in the Step 3 UI tab to skip OptiType entirely.

---

## Running locally (without Docker)

Requires Python 3.12+, Java 17+, and system packages: `bwa samtools bcftools tabix vienna-rna`.

```bash
# Install Python dependencies
pip install uv
uv sync

# Launch UI
.venv/bin/streamlit run app.py
```

---

## Output

Each run produces a versioned directory under `results/`:

```
results/run_<sample>_<timestamp>/
├── step1/   sorted.bam
├── step2/   filtered_variants.vcf.gz
├── step3/   hla_alleles.txt
├── step4/   candidate_neoantigens.tsv
├── step5/   ranked_neoantigens.tsv
├── step6/   ordered_epitopes.fasta
├── step7/   vaccine_mrna_*.fasta  candidate_comparison.json
└── step9/   vaccine_report.md  figures/
```

---

## Status

- [x] Step 1: Preprocessing
- [x] Step 2: Variant Calling
- [x] Step 3: HLA Typing
- [x] Step 4: Neoantigen Prediction
- [x] Step 5: Immunogenicity Ranking
- [x] Step 6: Epitope Ordering
- [x] Step 7: mRNA Design
- [x] Step 9: Report Generation
- [ ] Snakemake workflow (headless, no UI)
- [ ] Per-step Docker containers
- [ ] Full biological end-to-end validation on real tumor/normal data

---

## Disclaimer

This is a research proof of concept. It is not validated for clinical use and should not be used to guide medical decisions.
