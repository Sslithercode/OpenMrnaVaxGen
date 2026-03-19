# mRNA Vaccine Pipeline

> **Research use only. Not for clinical or diagnostic use.**

An open-source computational pipeline for designing personalized mRNA cancer vaccines. Given tumor/normal sequencing data and a patient's HLA profile, it identifies tumor-specific neoantigens, ranks them by predicted immunogenicity, and encodes the top candidates into an optimized mRNA vaccine sequence ready for wet lab synthesis.

## Why this exists

Every cancer patient's tumor carries a unique set of somatic mutations. Some of those mutations produce altered proteins that the immune system can, in principle, recognize as foreign — these are called **neoantigens**. Which neoantigens are visible to the immune system depends on the patient's HLA alleles, which determine what peptide fragments get displayed on the cell surface.

This pipeline automates the full computational workflow: calling somatic variants from tumor/normal sequencing, predicting which mutant peptides bind the patient's specific HLA alleles, ranking candidates by predicted immunogenicity, and assembling the top epitopes into an optimized mRNA construct. The output is a FASTA file and synthesis report you can hand to a molecular biology core facility.

### What a real run produces

Running on HCC1143 (a triple-negative breast cancer cell line), chr17 only, the pipeline:

- Called somatic variants with GATK Mutect2 (tumor vs. matched normal)
- Predicted MHC binding for all candidate peptides with MHCflurry 2.0
- Ranked 30 neoantigen candidates from ~158 passing filters (128 removed as subclonal, VAF < 0.05)
- Top candidate: **RSVAGVLHR** — HLA-A\*31:01, 25.7 nM affinity, presentation score 0.907, agretopicity 2.44
- Assembled 13 unique epitopes into a 179 aa polyepitope string joined with GPGPG linkers
- Codon-optimized and assembled a **870 nt mRNA construct** (GC 52.9%, MFE −348 kcal/mol)

---

## Pipeline

```
Input: Tumor WES + Normal WES
       ↓
Step 1: Preprocessing          → Clean BAMs             [GATK MarkDuplicates + SortSam]
Step 2: Variant Calling        → Somatic VCF            [GATK Mutect2 + FilterMutectCalls]
Step 3: HLA Typing             → HLA alleles            [OptiType]
Step 4: Neoantigen Prediction  → Candidate peptides     [MHCflurry 2.0]
Step 5: Immunogenicity Ranking → Ranked neoantigen list [composite score, IMPROVE weights]
Step 6: Epitope Ordering       → Ordered epitope string [Held-Karp exact TSP / greedy]
Step 7: mRNA Design            → Full mRNA sequence     [VaxPress + LinearFold/RNAfold]
Step 9: Report                 → vaccine_report.md
       ↓
Output: results/run_<sample>_<timestamp>/
```

**Entry A — full:** Raw FASTQs / BAMs → Steps 1 → 2 → 3 → 4 → 5 → 6 → 7 → 9

**Entry B — pre-called:** Existing MAF/VCF → Steps 3 → 4 → 5 → 6 → 7 → 9

---

## Prerequisites

The following are not versioned and must be set up manually.

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
| `hg38.fa` + `hg38.fa.fai` + `hg38.dict` | Steps 2, 4 (Mutect2 / MHCflurry) |
| `b37.20.21.fasta` + index | Step 1 (preprocessing) |

The hg38 reference is available from GATK's Google Cloud bucket (`gs://genomics-public-data/resources/broad/hg38/`) or UCSC. Both require indexing with `samtools faidx` and `gatk CreateSequenceDictionary`.

### 4. Add input data

Drop your tumor and normal BAM/BAI files into `data/test/`:

```
data/test/
├── test_tumor.bam
├── test_tumor.bam.bai
├── normal_chr17.bam
└── normal_chr17.bam.bai
```

---

## Running with Docker (experimental)

### Launch the Streamlit UI

```bash
docker compose up app
```

Open [http://localhost:8501](http://localhost:8501). Configure each step in its tab and run them in order.

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

### Step 3: HLA Typing (FALLBACK)


### Preferred
If you already have HLA alleles typed from clinical sequencing, use the **Manual HLA entry** panel in the Step 3 UI tab to skip OptiType entirely.


HLA typing uses OptiType, which runs as a separate container between two phases.

**Phase 1** — extract HLA reads:
```bash
docker compose run --rm pipeline python3 scripts/hla_typing.py --extract
```

**Phase 2** — run OptiType:
```bash
OPTITYPE_DATA_DIR=./results/run_<id>/step3 \
OPTITYPE_SAMPLE=hcc1143_normal \
docker compose run --rm optitype
```

**Phase 3** — parse results:
```bash
docker compose run --rm pipeline python3 scripts/hla_typing.py --parse
```


---

## Running locally (without Docker)

Requires Python 3.12+, Java 17+, and system packages: `bwa samtools bcftools tabix vienna-rna`.

```bash
pip install uv
uv sync
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
├── step5/   ranked_neoantigens.tsv  all_scored_candidates.tsv  subclonal_filtered.tsv
├── step6/   ordered_epitopes.fasta  junction_scores.tsv
├── step7/   vaccine_mrna_*.fasta  candidate_comparison.json
└── step9/   vaccine_report.md  figures/
```

The final report (`vaccine_report.md`) includes ranked candidate tables, HLA coverage plots, junction score heatmaps, mRNA construct metrics, and a wet lab synthesis checklist.

---

## Ranking methodology

Candidates are scored using a composite immunogenicity metric based on feature importance from [IMPROVE (Frontiers Immunology, 2024)](https://doi.org/10.3389/fimmu.2024.1397590):

these parameters are adjustable in the streamlit ui.

| Feature | Weight | Rationale |
|---------|--------|-----------|
| MHCflurry presentation score | 40% | Most predictive of actual antigen presentation |
| Agretopicity (DAI) | 25% | log2(WT affinity / mutant affinity) — higher means T cell repertoire not tolerized to this peptide |
| VAF / clonality | 20% | Clonal mutations present in more tumor cells; clinical standard requires VAF ≥ 0.05 |
| BLOSUM mutation score | 10% | More radical amino acid substitution = more foreign to immune system |
| Foreignness | 5% | BLOSUM kernel similarity between mutant and wildtype peptide |

Candidates are flagged (but not removed) for `NEG_AGRETOPICITY` (mutant binds MHC worse than wildtype) and `LOW_VAF` (0.05–0.10, moderately subclonal).

---

## Epitope ordering

To minimize immunogenic junctional peptides at GPGPG linker boundaries, epitope concatenation order is solved as a Hamiltonian path problem where edge weights are the worst MHCflurry presentation score across all junction peptides between each pair of epitopes. The pipeline uses:

- **Held-Karp exact dynamic programming** for N ≤ 15 epitopes — globally optimal
- **Greedy nearest-neighbor** for N > 15 — O(N²), good approximation

---

## Known limitations

- **Chromosome scope**: The default test data covers chr17 only. Genome-wide variant calling will produce substantially more candidates.
- **Cell line vs. primary tumor**: HCC1143 is a cell line. Clonal architecture and HLA expression differ from patient tumors.
- **No TCR validation**: Candidates are selected on MHC binding/presentation scores only. NetTCR-2.2 or IEDB T cell immunogenicity scoring is not currently integrated.
- **MFE/nt metric**: The LinearDesign optimal range (−0.48 to −0.60 kcal/mol/nt) was derived from full-length protein antigens; short polyepitope constructs (~870 nt) will fall outside this range due to UTR dilution — this is expected, not a failure.
- **CodonFM validation**: Step 8 (nvidia/NV-CodonFM-80M scoring) is implemented as a placeholder. Empirical HEK293T expression testing is the recommended ranking method for now.

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
- [ ] Step 8: CodonFM validation
- [ ] Snakemake headless workflow
- [ ] Per-step Docker containers
- [ ] Genome-wide variant calling (currently chr17 toest data)
- [ ] NetTCR-2.2 / IEDB T cell immunogenicity integration

---

### RoadMap
- [ ] Run benchmarks starting from neantigen generation step with public data
- [ ] Run 1 genome wide test
- [ ] NetTCR-2.2 / IEDB T cell immunogenicity integration

## License

MIT

LinearFold is provided under a Non Commercial license: if you want this fully open source use vienna
---

## Citation

If you use this pipeline in your research, please cite the tools it depends on:

- **MHCflurry**: O'Donnell et al., *Cell Systems*, 2020
- **VaxPress**: Ju, Ku & Chang, 2023
- **GATK/Mutect2**: Van der Auwera & O'Connor, *Bioinformatics Data Skills*, 2020
- **OptiType**: Szolek et al., *Bioinformatics*, 2014
- **LinearFold**: Huang et al., *ISMB*, 2019
- **IMPROVE ranking weights**: Frontiers in Immunology, 2024
