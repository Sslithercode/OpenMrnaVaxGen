# Melanoma mRNA Vaccine Pipeline — Project Plan

## Goal
Build an open-source computational pipeline to design personalized mRNA vaccines for melanoma.
Given tumor/normal sequencing data and a patient's HLA profile, identify tumor-specific neoantigens,
rank them by immunogenicity, and encode the top candidates into an optimized mRNA vaccine sequence.

---

## Pipeline Overview

```
Input: Tumor WES + Normal WES + (optional) RNA-seq
         ↓
Step 1: Preprocessing          → Clean BAMs
Step 2: Variant Calling        → Somatic VCF
Step 3: HLA Typing             → HLA alleles
Step 4: Neoantigen Prediction  → Candidate peptides + MHC binding scores
Step 5: Immunogenicity Ranking → Ranked neoantigen list
Step 6: Epitope Ordering       → Ordered epitope string
Step 7: mRNA Design            → Full mRNA sequence (FASTA)
Step 8: Validation & Report    → Quality metrics + PDF report

Output: vaccine_mrna.fasta + vaccine_design_report.pdf
```

---

## Entry Points

The pipeline supports two entry points:

- **Entry A (full):** Raw FASTQs → Steps 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8
- **Entry B (pre-called):** MAF/VCF file → Steps 3 → 4 → 5 → 6 → 7 → 8

---

## Steps Detail

### Step 1: Preprocessing
**Goal:** Convert raw sequencing reads into clean, analysis-ready BAMs

**Tools:** BWA-MEM, GATK MarkDuplicates, GATK SortSam
**Inputs:** Tumor FASTQ (R1, R2), Normal FASTQ (R1, R2), Reference genome (hg38)
**Outputs:** tumor.final.bam, normal.final.bam
**Status:** ✅ Script written and validated on test data (NA12878 chr20/21)

---

### Step 2: Somatic Variant Calling
**Goal:** Identify mutations present in tumor but not normal (the "diff")

**Tools:** GATK Mutect2, GATK FilterMutectCalls, VEP (annotation)
**Inputs:** tumor.final.bam, normal.final.bam, reference genome
**Outputs:** filtered_variants.vcf.gz
**Status:** ✅ Script written and validated. Note: real tumor/normal test data needed for
full biological validation. Current test data (NA12878 synthetic normal) returns 0 variants
as expected — pipeline code is correct.

**Pending:** Download HCC1143 breast cancer tumor/normal BAMs (real cancer data) for proper
end-to-end validation. Files are on Google Cloud Storage (gs://gatk-best-practices/somatic-hg38/).
Requires gsutil + Google account.

---

### Step 3: HLA Typing
**Goal:** Determine patient's HLA alleles — the "display protocol" of their immune system

**Tools:** OptiType
**Inputs:** normal.final.bam
**Outputs:** HLA allele calls (e.g. HLA-A*02:01, HLA-B*07:02, HLA-C*05:01)
**Status:** ✅ Script written (scripts/hla_typing.py)

---

### Step 4: Neoantigen Prediction & MHC Binding
**Goal:** Generate candidate peptides from mutations and predict which bind the patient's HLA

**Tools:** pVACseq, MHCflurry
**Inputs:** filtered_variants.vcf.gz, HLA alleles
**Outputs:** candidate_neoantigens.tsv (peptide, HLA, IC50 binding score)
**Status:** 🔲 Not started

---

### Step 5: Immunogenicity Scoring & Ranking
**Goal:** Score candidates beyond binding affinity — predict actual T-cell response

**Tools:** BigMHC, AlphaFold (optional, structural divergence score)
**Inputs:** candidate_neoantigens.tsv, RNA-seq TPM values (optional), VAF from VCF
**Scoring formula:**
```
score = w1 * (1/ic50) + w2 * immunogenicity + w3 * expression_tpm + w4 * vaf
```
**Outputs:** ranked_neoantigens.tsv (top 20-34 candidates)
**Status:** 🔲 Not started

---

### Step 6: Epitope Ordering
**Goal:** Find optimal concatenation order to minimize junctional neoantigens

**Tools:** pVACvector (graph optimization / TSP variant)
**Inputs:** ranked_neoantigens.tsv
**Outputs:** ordered_epitope_string.fasta (epitopes joined with GPGPG linkers)
**Status:** 🔲 Not started

---

### Step 7: mRNA Sequence Design
**Goal:** Design a full optimized mRNA construct encoding the epitope string

**Tools:** LinearDesign (codon optimization + secondary structure), Optimus 5-Prime (UTR scoring), RNAfold (validation)
**Components of final mRNA:**
- 5' cap (m7G)
- 5' UTR (selected from validated library via Optimus 5-Prime)
- Kozak sequence
- Start codon (ATG)
- Epitope coding sequence (LinearDesign optimized)
- 3' UTR (beta-globin standard)
- Poly-A tail (120 nt)
- m1Ψ substitution recommendation (wet lab flag)

**Outputs:** vaccine_mrna.fasta
**Status:** 🔲 Not started

---

### Step 8: Validation & Report
**Goal:** Validate mRNA quality and generate human-readable report

**Tools:** RNAfold (MFE score), custom Python, ReportLab (PDF)
**Checks:**
- Minimum free energy (MFE) — lower = more stable secondary structure
- No accidental ORFs
- No unwanted splice sites
- No restriction enzyme sites

**Outputs:** vaccine_design_report.pdf, quality_metrics.json
**Status:** 🔲 Not started

---

## Test Data

| Dataset | Purpose | Status |
|---|---|---|
| NA12878 chr20/21 BAM | Steps 1-2 code validation | ✅ Downloaded |
| Synthetic normal (wgsim) | Steps 1-2 code validation | ✅ Generated |
| TCGA-SKCM MAF (469 patients) | Steps 4-8 real melanoma data | ✅ Downloaded |
| HCC1143 tumor/normal BAMs | Steps 1-2 biological validation | 🔲 Pending (gsutil) |

---

## Tools Stack

| Step | Tool | License | Install |
|---|---|---|---|
| Alignment | BWA-MEM | MIT | apt |
| BAM processing | GATK/Picard | BSD | manual |
| Variant calling | GATK Mutect2 | BSD | manual |
| HLA typing | OptiType | MIT | uv |
| MHC binding | MHCflurry | Apache 2.0 | uv |
| Immunogenicity | BigMHC | MIT | uv |
| Protein structure | AlphaFold (optional) | Apache 2.0 | separate |
| Epitope ordering | pVACtools | BSD | uv |
| Codon optimization | LinearDesign | Academic | manual |
| UTR scoring | Optimus 5-Prime | MIT | uv |
| RNA structure | RNAfold (ViennaRNA) | MIT | apt |
| Reporting | ReportLab | BSD | uv |

---

## Open Source Strategy

- GitHub repo: `melanoma-vaccine-pipeline`
- Workflow manager: Snakemake (ties all steps together)
- Containerization: Docker (one container per step for reproducibility)
- Documentation: ReadTheDocs
- Test data: TCGA-SKCM MAF (open access) + synthetic test BAMs

---

## Current Status

- [x] Project structure created
- [x] System dependencies installed (BWA, samtools, GATK, bcftools)
- [x] Python environment set up (uv + pyproject.toml)
- [x] Step 1 script written and validated
- [x] Step 2 script written and validated
- [x] TCGA-SKCM MAF downloaded (real melanoma data for Steps 4-8)
- [ ] HCC1143 real tumor/normal BAMs (biological validation of Steps 1-2)
- [x] Step 3: HLA typing
- [ ] Step 4: Neoantigen prediction
- [ ] Step 5: Immunogenicity ranking
- [ ] Step 6: Epitope ordering
- [ ] Step 7: mRNA design
- [ ] Step 8: Validation and report
- [ ] Snakemake workflow tying all steps together
- [ ] Docker containers
- [ ] GitHub repo + documentation
