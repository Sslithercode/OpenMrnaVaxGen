# Melanoma mRNA Vaccine Pipeline

> **Work in progress / proof of concept.** Not for clinical use.

An open-source computational pipeline for designing personalized mRNA vaccines for melanoma. Given tumor/normal sequencing data and a patient's HLA profile, the pipeline identifies tumor-specific neoantigens, ranks them by immunogenicity, and encodes the top candidates into an optimized mRNA vaccine sequence.

The UI is a [Streamlit](https://streamlit.io) app that lets you run each step interactively with configurable parameters.

---

## Pipeline Overview

```
Input: Tumor WES + Normal WES + (optional) RNA-seq
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

**Entry A (full):** Raw FASTQs → Steps 1 → 2 → 3 → 4 → 5 → 6 → 7 → 9

**Entry B (pre-called):** MAF/VCF → Steps 3 → 4 → 5 → 6 → 7 → 9

---

## Tools

| Step | Tool | License |
|------|------|---------|
| BAM processing | GATK/Picard | BSD |
| Variant calling | GATK Mutect2 | BSD |
| HLA typing | OptiType | MIT |
| MHC binding | MHCflurry | Apache 2.0 |
| Epitope ordering | TSP / greedy (custom) | — |
| Codon optimization | VaxPress | MIT |
| RNA structure | RNAfold (ViennaRNA) | MIT |

---

## Requirements

- Python 3.12+
- Java (for GATK)
- [uv](https://github.com/astral-sh/uv) (package manager)
- BWA, samtools, bcftools (system packages)
- GATK 4.5.0.0 (bundled under `tools/`)
- OptiType (Docker recommended for HLA typing step)
- Reference genomes: `reference/hg38.fa`, `reference/b37.20.21.fasta`

---

## Setup

```bash
git clone https://github.com/your-username/melanoma-pipeline
cd melanoma-pipeline

# Install Python dependencies
uv sync

# Launch the UI
.venv/bin/streamlit run app.py
```

---

## Usage

1. Open the Streamlit UI in your browser.
2. Select an entry point in the sidebar (Entry A for raw FASTQs, Entry B if you already have a VCF/MAF).
3. Configure paths and parameters in each step's tab.
4. Run steps sequentially. Each run gets its own versioned output directory (`results/run_<sample>_<timestamp>/`).

For HLA typing (Step 3), OptiType must be run manually between Phase 1 (FASTQ extraction) and Phase 3 (result parsing). A manual HLA entry option is also available if you already have alleles.

---

## Output

Each run produces a directory under `results/` containing per-step subdirectories and a final `vaccine_report.md` with:

- Ranked neoantigen table
- mRNA construct sequence (FASTA)
- Candidate comparison (GC%, MFE, length)
- Quality metrics plots

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
- [ ] Snakemake workflow (command-line, no UI)
- [ ] Docker containers
- [ ] Biological end-to-end validation on real tumor/normal data

---

## Disclaimer

This is a research proof of concept. It is not validated for clinical use and should not be used to guide medical decisions.
