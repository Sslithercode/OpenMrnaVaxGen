"""
Melanoma mRNA Vaccine Pipeline — Streamlit UI
Orchestrates all pipeline steps with adjustable parameters.
"""

import os
import json
import sys
import tempfile
import textwrap
import subprocess
from pathlib import Path

import streamlit as st

# ── Paths / run resolver ───────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).parent
VENV_PYTHON  = PROJECT_ROOT / ".venv" / "bin" / "python"
SCRIPTS_DIR  = PROJECT_ROOT / "scripts"
TOOLS_DIR    = PROJECT_ROOT / "tools"
REFERENCE_DIR = PROJECT_ROOT / "reference"
DATA_DIR     = PROJECT_ROOT / "data"

sys.path.insert(0, str(SCRIPTS_DIR))
import scripts.paths as _paths  # noqa: E402

STEP_NAMES = {
    1: "Preprocessing",
    2: "Variant Calling",
    3: "HLA Typing",
    4: "Neoantigen Prediction",
    5: "Immunogenicity Ranking",
    6: "Epitope Ordering",
    7: "mRNA Design",
    9: "Report",
}

ENTRY_B_STEPS = {3, 4, 5, 6, 7, 9}  # steps available when starting from VCF

# ── Session state init ────────────────────────────────────────────────────────

ALL_STEPS = [1, 2, 3, 4, 5, 6, 7, 9]

def init_state():
    if "run_id" not in st.session_state:
        st.session_state.run_id = _paths.RUN_ID
    if "step_status" not in st.session_state:
        st.session_state.step_status = {i: "pending" for i in ALL_STEPS}
    if "step_logs" not in st.session_state:
        st.session_state.step_logs = {i: "" for i in ALL_STEPS}


def active_run_dir() -> Path:
    return _paths.RESULTS / st.session_state.run_id


def step_default(n: int) -> Path:
    """Default output directory for step n in the active run."""
    return active_run_dir() / f"step{n}"

# ── Runner ────────────────────────────────────────────────────────────────────

def run_script(script_code: str, log_placeholder, step_num: int) -> bool:
    """Write a Python script to a temp file, run it with the venv Python,
    and stream stdout+stderr to the given Streamlit placeholder.
    Passes PIPELINE_RUN so all steps write into the active run directory."""
    with tempfile.NamedTemporaryFile(suffix=".py", mode="w", delete=False,
                                     prefix="pipeline_step_") as f:
        f.write(script_code)
        tmp_path = f.name

    st.session_state.step_status[step_num] = "running"
    output_lines = []

    env = os.environ.copy()
    env["PIPELINE_RUN"] = st.session_state.run_id

    try:
        process = subprocess.Popen(
            [str(VENV_PYTHON), tmp_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=str(PROJECT_ROOT),
            env=env,
        )

        for line in process.stdout:
            output_lines.append(line.rstrip())
            log_placeholder.code("\n".join(output_lines[-200:]), language="")

        process.wait()
        st.session_state.step_logs[step_num] = "\n".join(output_lines)

        success = process.returncode == 0
        st.session_state.step_status[step_num] = "done" if success else "failed"
        return success

    except Exception as e:
        st.session_state.step_status[step_num] = "failed"
        st.session_state.step_logs[step_num] = str(e)
        log_placeholder.code(str(e), language="")
        return False

    finally:
        os.unlink(tmp_path)


def step_header(step_num: int):
    status = st.session_state.step_status[step_num]
    icons = {"pending": "⬜", "running": "🔄", "done": "✅", "failed": "❌"}
    st.markdown(f"### {icons[status]} Step {step_num}: {STEP_NAMES[step_num]}")


def show_prior_log(step_num: int):
    log = st.session_state.step_logs.get(step_num, "")
    if log:
        with st.expander("Previous run output", expanded=False):
            st.code(log, language="")

# ── Page config ───────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="Melanoma mRNA Vaccine Pipeline",
    layout="wide",
    initial_sidebar_state="expanded",
)

init_state()

# ── Sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    st.title("🧬 Pipeline Config")

    entry_point = st.radio(
        "Entry Point",
        ["A — Full (FASTQs → mRNA)", "B — Pre-called (VCF/MAF → mRNA)"],
        key="entry_point",
    )
    entry_b = entry_point.startswith("B")

    st.divider()
    st.subheader("Active Run")
    st.code(st.session_state.run_id, language="")
    st.caption(str(active_run_dir()))

    if st.button("🆕 New Run", help="Start a fresh versioned output directory"):
        new_id = _paths.new_run()
        st.session_state.run_id = new_id
        st.session_state.step_status = {i: "pending" for i in ALL_STEPS}
        st.session_state.step_logs   = {i: "" for i in ALL_STEPS}
        st.rerun()

    past_runs = sorted(
        [d.name for d in _paths.RESULTS.iterdir() if d.is_dir() and d.name.startswith("run_")],
        reverse=True,
    ) if _paths.RESULTS.exists() else []
    if past_runs:
        selected = st.selectbox("Switch to past run", ["(current)"] + past_runs, key="past_run_select")
        if selected != "(current)" and selected != st.session_state.run_id:
            if st.button("Load selected run"):
                st.session_state.run_id = selected
                st.session_state.step_status = {i: "pending" for i in ALL_STEPS}
                st.session_state.step_logs   = {i: "" for i in ALL_STEPS}
                st.rerun()

    st.divider()
    st.subheader("Shared Paths")

    gatk_jar = st.text_input(
        "GATK JAR",
        str(TOOLS_DIR / "gatk-4.5.0.0" / "gatk-package-4.5.0.0-local.jar"),
        key="gatk_jar",
    )
    ref_hg38 = st.text_input(
        "Reference genome (hg38)",
        str(REFERENCE_DIR / "hg38.fa"),
        key="ref_hg38",
    )
    ref_b37 = st.text_input(
        "Reference genome (b37)",
        str(REFERENCE_DIR / "b37.20.21.fasta"),
        key="ref_b37",
    )

    st.divider()
    st.subheader("Step Status")
    icons = {"pending": "⬜", "running": "🔄", "done": "✅", "failed": "❌"}
    for i in ALL_STEPS:
        if entry_b and i not in ENTRY_B_STEPS:
            continue
        status = st.session_state.step_status[i]
        st.write(f"{icons[status]} Step {i}: {STEP_NAMES[i]}")

    if st.button("Reset all statuses"):
        st.session_state.step_status = {i: "pending" for i in ALL_STEPS}
        st.session_state.step_logs   = {i: "" for i in ALL_STEPS}
        st.rerun()

# ── Main tabs ─────────────────────────────────────────────────────────────────

tab_labels = (
    ["Overview"]
    + ([] if entry_b else ["1: Preprocess", "2: Variants"])
    + ["3: HLA", "4: Neoantigens", "5: Rank", "6: Order", "7: mRNA", "9: Report"]
)
tabs = st.tabs(tab_labels)
tab_idx = 0

# ─── Overview ─────────────────────────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1

    st.header("Melanoma mRNA Vaccine Pipeline")
    st.markdown(
        """
        This pipeline takes tumor/normal sequencing data and a patient HLA profile,
        identifies tumor-specific neoantigens, ranks them by immunogenicity, and encodes
        the top candidates into an optimized mRNA vaccine sequence.

        **Entry A (full):** Raw FASTQs → Steps 1 → 2 → 3 → 4 → 5 → 6 → 7 → 9

        **Entry B (pre-called):** MAF/VCF → Steps 3 → 4 → 5 → 6 → 7 → 9

        Each run generates its own versioned output directory. Select the entry point
        in the sidebar and configure per-step parameters in each tab.
        """
    )

    st.info(f"**Active run:** `{st.session_state.run_id}`  \n**Output:** `{active_run_dir()}`")

    st.subheader("Pipeline Diagram")
    st.code(
        """
Input: Tumor WES + Normal WES + (optional) RNA-seq
       ↓
Step 1: Preprocessing          → Clean BAMs             [GATK MarkDuplicates + SortSam]
Step 2: Variant Calling        → Somatic VCF            [GATK Mutect2 + FilterMutectCalls]
Step 3: HLA Typing             → HLA alleles            [OptiType]
Step 4: Neoantigen Prediction  → Candidate peptides     [MHCflurry]
Step 5: Immunogenicity Ranking → Ranked neoantigen list [composite score]
Step 6: Epitope Ordering       → Ordered epitope string [TSP/greedy]
Step 7: mRNA Design            → Full mRNA sequence     [VaxPress + RNAfold]
Step 9: Report                 → vaccine_report.md      [matplotlib + markdown]
       ↓
Output: results/run_<sample>_<timestamp>/
        """,
        language="",
    )

    st.subheader("Current Status")
    cols = st.columns(len(ALL_STEPS))
    icons = {"pending": "⬜", "running": "🔄", "done": "✅", "failed": "❌"}
    for idx, i in enumerate(ALL_STEPS):
        with cols[idx]:
            status = st.session_state.step_status[i]
            st.metric(f"Step {i}", icons[status])

# ─── Step 1: Preprocessing ────────────────────────────────────────────────────

if not entry_b:
    with tabs[tab_idx]:
        tab_idx += 1
        step_header(1)
        st.markdown(
            "Convert raw sequencing reads into clean, analysis-ready BAMs using "
            "GATK MarkDuplicates and SortSam."
        )

        with st.expander("Parameters", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                s1_input_bam = st.text_input(
                    "Input BAM",
                    str(DATA_DIR / "test" / "test_tumor.bam"),
                    key="s1_input_bam",
                )
                s1_out_dir = st.text_input(
                    "Output directory",
                    str(step_default(1)),
                    key="s1_out_dir",
                )
            with col2:
                s1_gatk_jar = st.text_input("GATK JAR", gatk_jar, key="s1_gatk_jar")
                s1_reference = st.text_input(
                    "Reference", ref_b37, key="s1_reference"
                )

        show_prior_log(1)
        log1 = st.empty()

        if st.button("Run Step 1: Preprocessing", type="primary"):
            script = textwrap.dedent(f"""\
                import sys
                sys.path.insert(0, {str(SCRIPTS_DIR)!r})
                import preprocess
                from pathlib import Path
                preprocess.GATK_JAR  = Path({s1_gatk_jar!r})
                preprocess.REFERENCE = Path({s1_reference!r})
                preprocess.INPUT_BAM = Path({s1_input_bam!r})
                preprocess.OUT_DIR   = Path({s1_out_dir!r})
                preprocess.main()
            """)
            run_script(script, log1, 1)
            st.rerun()

# ─── Step 2: Variant Calling ──────────────────────────────────────────────────

if not entry_b:
    with tabs[tab_idx]:
        tab_idx += 1
        step_header(2)
        st.markdown(
            "Identify somatic mutations present in tumor but not normal using "
            "GATK Mutect2 and FilterMutectCalls."
        )

        with st.expander("Parameters", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                s2_tumor_bam = st.text_input(
                    "Tumor BAM",
                    str(step_default(1) / "sorted.bam"),
                    key="s2_tumor_bam",
                )
                s2_normal_bam = st.text_input(
                    "Normal BAM",
                    str(DATA_DIR / "test" / "normal_chr17.bam"),
                    key="s2_normal_bam",
                )
                s2_out_dir = st.text_input(
                    "Output directory",
                    str(step_default(2)),
                    key="s2_out_dir",
                )
            with col2:
                s2_gatk_jar  = st.text_input("GATK JAR", gatk_jar, key="s2_gatk_jar")
                s2_reference = st.text_input("Reference", ref_hg38, key="s2_reference")
                s2_tumor_name  = st.text_input("Tumor sample name",  "HCC1143_tumor",  key="s2_tumor_name")
                s2_normal_name = st.text_input("Normal sample name", "HCC1143_normal", key="s2_normal_name")
                s2_intervals   = st.text_input("Intervals (chromosome)", "chr17", key="s2_intervals")

        show_prior_log(2)
        log2 = st.empty()

        if st.button("Run Step 2: Variant Calling", type="primary"):
            script = textwrap.dedent(f"""\
                import sys
                sys.path.insert(0, {str(SCRIPTS_DIR)!r})
                import variant
                from pathlib import Path
                variant.GATK_JAR   = Path({s2_gatk_jar!r})
                variant.REFERENCE  = Path({s2_reference!r})
                variant.TUMOR_BAM  = Path({s2_tumor_bam!r})
                variant.NORMAL_BAM = Path({s2_normal_bam!r})
                variant.OUT_DIR    = Path({s2_out_dir!r})
                # patch sample names and intervals inside the Mutect2 command
                _orig_main = variant.main
                def _patched_main():
                    import subprocess, sys as _sys
                    out = variant.OUT_DIR
                    out.mkdir(parents=True, exist_ok=True)
                    raw_vcf      = out / "raw_variants.vcf.gz"
                    filtered_vcf = out / "filtered_variants.vcf.gz"
                    variant.run(variant.gatk([
                        "Mutect2",
                        "-R",          str(variant.REFERENCE),
                        "-I",          str(variant.TUMOR_BAM),
                        "-I",          str(variant.NORMAL_BAM),
                        "-tumor",      {s2_tumor_name!r},
                        "-normal",     {s2_normal_name!r},
                        "-O",          str(raw_vcf),
                        "--intervals", {s2_intervals!r},
                    ]), "Mutect2")
                    variant.run(variant.gatk([
                        "FilterMutectCalls",
                        "-R", str(variant.REFERENCE),
                        "-V", str(raw_vcf),
                        "-O", str(filtered_vcf),
                    ]), "FilterMutectCalls")
                    result = subprocess.run(
                        ["java", "-jar", str(variant.GATK_JAR), "CountVariants", "-V", str(filtered_vcf)],
                        capture_output=True, text=True
                    )
                    print("── Variant Count ──")
                    print(result.stdout)
                    print(f"Step 2 complete. Output: {{filtered_vcf}}")
                variant.main = _patched_main
                variant.main()
            """)
            run_script(script, log2, 2)
            st.rerun()

# ─── Step 3: HLA Typing ───────────────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(3)
    st.markdown(
        "Determine the patient's HLA alleles using OptiType. "
        "Runs in two phases: extract reads, then parse OptiType results. "
        "OptiType must be run between the two phases (see sidebar for Docker command)."
    )

    with st.expander("Manual HLA entry (skip OptiType)", expanded=False):
        st.markdown("If you already have HLA alleles, enter them below to write `hla_alleles.txt` directly.")
        manual_hla = st.text_area(
            "HLA alleles (one per line, e.g. HLA-A*02:01)",
            "HLA-A*02:01\nHLA-A*24:02\nHLA-B*07:02\nHLA-B*35:01\nHLA-C*04:01\nHLA-C*07:02",
            key="manual_hla",
        )
        s3_out_dir_manual = st.text_input(
            "Output directory (for hla_alleles.txt)",
            str(step_default(3)),
            key="s3_out_dir_manual",
        )
        if st.button("Write HLA alleles manually"):
            out = Path(s3_out_dir_manual)
            out.mkdir(parents=True, exist_ok=True)
            alleles = [a.strip() for a in manual_hla.strip().splitlines() if a.strip()]
            with open(out / "hla_alleles.txt", "w") as f:
                f.write("\n".join(alleles) + "\n")
            st.success(f"Wrote {len(alleles)} alleles to {out / 'hla_alleles.txt'}")
            st.session_state.step_status[3] = "done"

    st.divider()
    with st.expander("Phase 1 & 3 Parameters", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            s3_normal_bam = st.text_input(
                "Normal BAM",
                str(DATA_DIR / "test" / "normal_chr17.bam"),
                key="s3_normal_bam",
            )
            s3_out_dir = st.text_input(
                "Output directory",
                str(step_default(3)),
                key="s3_out_dir",
            )
        with col2:
            s3_sample_prefix = st.text_input("Sample prefix", "hcc1143_normal", key="s3_sample_prefix")
            s3_hla_region    = st.text_input(
                "HLA region (leave blank for all reads)",
                "",
                key="s3_hla_region",
                help="e.g. chr6:28510120-33480577 for full genome BAMs",
            )

    show_prior_log(3)
    log3 = st.empty()

    col_a, col_b = st.columns(2)
    with col_a:
        if st.button("Phase 1: Extract HLA reads", type="primary"):
            hla_region_val = s3_hla_region.strip() if s3_hla_region.strip() else "None"
            if hla_region_val != "None":
                hla_region_val = repr(s3_hla_region.strip())
            script = textwrap.dedent(f"""\
                import sys
                sys.path.insert(0, {str(SCRIPTS_DIR)!r})
                import hla_typing
                from pathlib import Path
                hla_typing.NORMAL_BAM    = Path({s3_normal_bam!r})
                hla_typing.OUT_DIR       = Path({s3_out_dir!r})
                hla_typing.SAMPLE_PREFIX = {s3_sample_prefix!r}
                hla_typing.HLA_REGION    = {hla_region_val}
                hla_typing.OUT_DIR.mkdir(parents=True, exist_ok=True)
                print("Step 3 — Phase 1: extracting HLA reads")
                hla_typing.extract_reads(hla_typing.NORMAL_BAM, hla_typing.OUT_DIR, hla_typing.HLA_REGION)
                print(f"Done. FASTQs written to {{hla_typing.OUT_DIR}}")
            """)
            run_script(script, log3, 3)
            st.rerun()

    with col_b:
        if st.button("Phase 3: Parse OptiType results"):
            script = textwrap.dedent(f"""\
                import sys
                sys.path.insert(0, {str(SCRIPTS_DIR)!r})
                import hla_typing
                from pathlib import Path
                hla_typing.NORMAL_BAM    = Path({s3_normal_bam!r})
                hla_typing.OUT_DIR       = Path({s3_out_dir!r})
                hla_typing.SAMPLE_PREFIX = {s3_sample_prefix!r}
                hla_typing.OUT_DIR.mkdir(parents=True, exist_ok=True)
                print("Step 3 — Phase 3: parsing OptiType results")
                hla_alleles = hla_typing.parse_results(hla_typing.OUT_DIR, hla_typing.SAMPLE_PREFIX)
                alleles_file = hla_typing.OUT_DIR / "hla_alleles.txt"
                with open(alleles_file, "w") as f:
                    for allele in hla_alleles:
                        f.write(allele + "\\n")
                print("── HLA Typing Results ──")
                for allele in hla_alleles:
                    print(f"  {{allele}}")
                print(f"Step 3 complete. Alleles saved to: {{alleles_file}}")
            """)
            run_script(script, log3, 3)
            st.rerun()

    st.info(
        "Between Phase 1 and Phase 3, run OptiType on the extracted FASTQs. "
        "Expected output: `results/step3/<sample_prefix>_result.tsv`"
    )

# ─── Step 4: Neoantigen Prediction ────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(4)
    st.markdown(
        "Generate candidate peptides from somatic mutations and predict MHC binding "
        "affinity using MHCflurry."
    )

    with st.expander("Parameters", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            s4_vcf = st.text_input(
                "VCF / MAF file",
                str(step_default(2) / "filtered_variants.vcf.gz"),
                key="s4_vcf",
            )
            s4_hla_file = st.text_input(
                "HLA alleles file (Step 3 output)",
                str(step_default(3) / "hla_alleles.txt"),
                key="s4_hla_file",
            )
            s4_out_dir = st.text_input(
                "Output directory",
                str(step_default(4)),
                key="s4_out_dir",
            )
        with col2:
            s4_reference = st.text_input("Reference", ref_hg38, key="s4_reference")
            s4_peptide_lengths = st.multiselect(
                "Peptide lengths",
                options=[8, 9, 10, 11, 12],
                default=[8, 9, 10, 11],
                key="s4_peptide_lengths",
            )
            s4_binding_cutoff = st.number_input(
                "Binding cutoff (IC50 nM)",
                min_value=50,
                max_value=5000,
                value=500,
                step=50,
                key="s4_binding_cutoff",
                help="MHC binding threshold — standard is 500 nM",
            )

    show_prior_log(4)
    log4 = st.empty()

    if st.button("Run Step 4: Neoantigen Prediction", type="primary"):
        script = textwrap.dedent(f"""\
            import sys
            sys.path.insert(0, {str(SCRIPTS_DIR)!r})
            import neoantigen_prediction
            from pathlib import Path
            neoantigen_prediction.VCF_FILE        = Path({s4_vcf!r})
            neoantigen_prediction.HLA_FILE        = Path({s4_hla_file!r})
            neoantigen_prediction.REFERENCE       = Path({s4_reference!r})
            neoantigen_prediction.OUT_DIR         = Path({s4_out_dir!r})
            neoantigen_prediction.PEPTIDE_LENGTHS = {s4_peptide_lengths!r}
            neoantigen_prediction.BINDING_CUTOFF  = {s4_binding_cutoff}
            neoantigen_prediction.main()
        """)
        run_script(script, log4, 4)
        st.rerun()

# ─── Step 5: Immunogenicity Ranking ───────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(5)
    st.markdown(
        "Score neoantigen candidates using a composite immunogenicity model "
        "(IMPROVE 2024 feature importance). Outputs ranked list of top vaccine targets."
    )

    with st.expander("Parameters", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            s5_input = st.text_input(
                "Step 4 output (candidate_neoantigens.tsv)",
                str(step_default(4) / "candidate_neoantigens.tsv"),
                key="s5_input",
            )
            s5_hla_file = st.text_input(
                "HLA alleles file",
                str(step_default(3) / "hla_alleles.txt"),
                key="s5_hla_file",
            )
            s5_out_dir = st.text_input(
                "Output directory",
                str(step_default(5)),
                key="s5_out_dir",
            )
            s5_top_n = st.number_input(
                "Top N candidates",
                min_value=5,
                max_value=50,
                value=30,
                step=1,
                key="s5_top_n",
            )
        with col2:
            st.markdown("**Composite score weights** (must sum to ~1.0)")
            st.caption("Based on IMPROVE (Frontiers Immunology 2024) feature importance")
            s5_w_pres = st.slider("Presentation score weight", 0.0, 1.0, 0.40, 0.05, key="s5_w_pres")
            s5_w_agr  = st.slider("Agretopicity (DAI) weight", 0.0, 1.0, 0.25, 0.05, key="s5_w_agr")
            s5_w_vaf  = st.slider("VAF / clonality weight",    0.0, 1.0, 0.20, 0.05, key="s5_w_vaf")
            s5_w_bl   = st.slider("BLOSUM mutation weight",    0.0, 1.0, 0.10, 0.05, key="s5_w_bl")
            s5_w_for  = st.slider("Foreignness weight",        0.0, 1.0, 0.05, 0.05, key="s5_w_for")
            total_w = s5_w_pres + s5_w_agr + s5_w_vaf + s5_w_bl + s5_w_for
            st.caption(f"Total weight: {total_w:.2f} {'✅' if abs(total_w - 1.0) < 0.01 else '⚠️ should be 1.0'}")

            st.markdown("**Filters**")
            s5_max_aff  = st.number_input("Max affinity (nM IC50)", 50, 5000, 500, 50, key="s5_max_aff")
            s5_min_pres = st.number_input("Min presentation score", 0.0, 1.0, 0.1, 0.05, key="s5_min_pres")

    show_prior_log(5)
    log5 = st.empty()

    if st.button("Run Step 5: Immunogenicity Ranking", type="primary"):
        script = textwrap.dedent(f"""\
            import sys
            sys.path.insert(0, {str(SCRIPTS_DIR)!r})
            import candidate_ranking
            from pathlib import Path
            candidate_ranking.STEP4_OUTPUT     = Path({s5_input!r})
            candidate_ranking.HLA_FILE         = Path({s5_hla_file!r})
            candidate_ranking.OUT_DIR          = Path({s5_out_dir!r})
            candidate_ranking.W_PRESENTATION   = {s5_w_pres}
            candidate_ranking.W_AGRETOPICITY   = {s5_w_agr}
            candidate_ranking.W_VAF            = {s5_w_vaf}
            candidate_ranking.W_BLOSUM         = {s5_w_bl}
            candidate_ranking.W_FOREIGNNESS    = {s5_w_for}
            candidate_ranking.MAX_AFFINITY_NM  = {s5_max_aff}
            candidate_ranking.MIN_PRESENTATION = {s5_min_pres}
            candidate_ranking.TOP_N            = {s5_top_n}
            candidate_ranking.main()
        """)
        run_script(script, log5, 5)
        st.rerun()

# ─── Step 6: Epitope Ordering ─────────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(6)
    st.markdown(
        "Find the optimal concatenation order of epitopes to minimize junctional "
        "neoantigens. Uses exact TSP (Held-Karp) for small N, greedy nearest-neighbor "
        "for large N."
    )

    with st.expander("Parameters", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            s6_input = st.text_input(
                "Step 5 output (ranked_neoantigens.tsv)",
                str(step_default(5) / "ranked_neoantigens.tsv"),
                key="s6_input",
            )
            s6_hla_file = st.text_input(
                "HLA alleles file",
                str(step_default(3) / "hla_alleles.txt"),
                key="s6_hla_file",
            )
            s6_out_dir = st.text_input(
                "Output directory",
                str(step_default(6)),
                key="s6_out_dir",
            )
        with col2:
            s6_linker       = st.text_input("GPGPG linker sequence", "GPGPG", key="s6_linker")
            s6_junction_len = st.number_input(
                "Junction peptide length",
                min_value=7, max_value=11, value=9, step=1, key="s6_junction_len",
                help="Length of junctional peptides to check for MHC binding",
            )
            s6_max_greedy = st.number_input(
                "Max epitopes for greedy (vs exact TSP)",
                min_value=5, max_value=50, value=30, step=1, key="s6_max_greedy",
                help="Above this count use greedy, below use exact TSP (2^N)",
            )

    show_prior_log(6)
    log6 = st.empty()

    if st.button("Run Step 6: Epitope Ordering", type="primary"):
        script = textwrap.dedent(f"""\
            import sys
            sys.path.insert(0, {str(SCRIPTS_DIR)!r})
            import epitope_ordering
            from pathlib import Path
            epitope_ordering.STEP5_OUTPUT  = Path({s6_input!r})
            epitope_ordering.HLA_FILE      = Path({s6_hla_file!r})
            epitope_ordering.OUT_DIR       = Path({s6_out_dir!r})
            epitope_ordering.LINKER        = {s6_linker!r}
            epitope_ordering.JUNCTION_LEN  = {s6_junction_len}
            epitope_ordering.MAX_GREEDY    = {s6_max_greedy}
            epitope_ordering.main()
        """)
        run_script(script, log6, 6)
        st.rerun()

# ─── Step 7: mRNA Design ──────────────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(7)
    st.markdown(
        "Design a full optimized mRNA vaccine construct encoding the epitope string. "
        "Runs VaxPress codon optimization with multiple candidate profiles, "
        "assembles the full construct, and validates with RNAfold."
    )

    with st.expander("Construct Parameters", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            s7_input = st.text_input(
                "Step 6 FASTA (ordered_epitopes.fasta)",
                str(step_default(6) / "ordered_epitopes.fasta"),
                key="s7_input",
            )
            s7_out_dir = st.text_input(
                "Output directory",
                str(step_default(7)),
                key="s7_out_dir",
            )
            s7_poly_a_len = st.number_input(
                "Poly-A tail length (nt)",
                min_value=60, max_value=300, value=120, step=10, key="s7_poly_a_len",
            )
        with col2:
            s7_vaxpress_iters = st.number_input(
                "VaxPress iterations",
                min_value=50, max_value=1000, value=200, step=50, key="s7_vaxpress_iters",
            )
            s7_vaxpress_procs = st.number_input(
                "VaxPress processes",
                min_value=1, max_value=32, value=8, step=1, key="s7_vaxpress_procs",
            )
            s7_utr5  = st.text_input("5' UTR", "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC", key="s7_utr5")
            s7_kozak = st.text_input("Kozak sequence", "GCCACCATG", key="s7_kozak")

    with st.expander("Candidate Profiles", expanded=True):
        st.markdown(
            "Each profile runs VaxPress with different weight overrides, sweeping the "
            "MFE/CAI tradeoff (mirroring LinearDesign Nature 2023 λ sweep). "
            "Step 8 CodonFM scoring ranks them."
        )

        s7_candidates = []
        profile_defaults = [
            ("A_balanced",    "Balanced MFE + CAI (default)", {}),
            ("B_stability",   "Stability-biased (longer half-life)", {"mfe": 10, "cai": 2, "loop": 5}),
            ("C_expression",  "Expression-biased (higher protein yield)", {"mfe": 2, "cai": 10, "gc": 3}),
            ("D_low_uridine", "Uridine-minimized (reduced innate immune activation)", {"ucount": 10, "mfe": 5, "cai": 3}),
        ]

        for i, (pid, pdesc, pweights) in enumerate(profile_defaults):
            with st.container():
                cols = st.columns([2, 3, 4])
                with cols[0]:
                    enabled = st.checkbox(f"Enable", value=True, key=f"s7_cand_enabled_{i}")
                    cand_id = st.text_input("ID", pid, key=f"s7_cand_id_{i}")
                with cols[1]:
                    cand_desc = st.text_input("Description", pdesc, key=f"s7_cand_desc_{i}")
                with cols[2]:
                    weight_keys = ["mfe", "cai", "ucount", "loop", "gc"]
                    weight_vals = {}
                    wcols = st.columns(5)
                    for j, wk in enumerate(weight_keys):
                        with wcols[j]:
                            default_val = pweights.get(wk, 0)
                            v = st.number_input(
                                wk, min_value=0, max_value=20,
                                value=default_val, step=1,
                                key=f"s7_w_{i}_{wk}",
                                label_visibility="visible",
                            )
                            if v > 0:
                                weight_vals[wk] = v
                if enabled:
                    s7_candidates.append({
                        "id": cand_id,
                        "description": cand_desc,
                        "weights": weight_vals,
                    })
                st.divider()

    show_prior_log(7)
    log7 = st.empty()

    if st.button("Run Step 7: mRNA Design", type="primary"):
        candidates_repr = repr(s7_candidates)
        script = textwrap.dedent(f"""\
            import sys
            sys.path.insert(0, {str(SCRIPTS_DIR)!r})
            import mrna_design
            from pathlib import Path
            mrna_design.STEP6_FASTA          = Path({s7_input!r})
            mrna_design.OUT_DIR              = Path({s7_out_dir!r})
            mrna_design.UTR5                 = {s7_utr5!r}
            mrna_design.KOZAK                = {s7_kozak!r}
            mrna_design.POLY_A_LENGTH        = {s7_poly_a_len}
            mrna_design.VAXPRESS_ITERATIONS  = {s7_vaxpress_iters}
            mrna_design.VAXPRESS_PROCESSES   = {s7_vaxpress_procs}
            mrna_design.CANDIDATES           = {candidates_repr}
            mrna_design.main()
        """)
        run_script(script, log7, 7)
        st.rerun()

    # Show output files if step is done
    if st.session_state.step_status[7] == "done":
        st.subheader("Output Files")
        out_dir = Path(st.session_state.get("s7_out_dir", str(step_default(7))))
        comparison_file = out_dir / "candidate_comparison.json"
        if comparison_file.exists():
            with open(comparison_file) as f:
                comparison = json.load(f)
            import pandas as pd
            rows = []
            for m in comparison:
                rows.append({
                    "Candidate": m["candidate_id"],
                    "Description": m["description"],
                    "GC %": f"{m['gc_content']:.1%}" if m.get("gc_content") else "N/A",
                    "MFE (kcal/mol)": f"{m['mfe_kcal_mol']:.2f}" if m.get("mfe_kcal_mol") else "N/A",
                    "Length (nt)": m.get("total_length_nt", "N/A"),
                    "QC": "ISSUES" if m.get("issues") else ("warn" if m.get("warnings") else "pass"),
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True)

# ─── Step 9: Report Generation ────────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(9)
    st.markdown(
        "Generate a Markdown report with dynamic tables and matplotlib plots "
        "from all pipeline outputs in the active run."
    )

    with st.expander("Parameters", expanded=True):
        s9_out_dir = st.text_input(
            "Output directory",
            str(step_default(9)),
            key="s9_out_dir",
        )

    show_prior_log(9)
    log9 = st.empty()

    if st.button("Run Step 9: Generate Report", type="primary"):
        script = textwrap.dedent(f"""\
            import sys
            sys.path.insert(0, {str(SCRIPTS_DIR)!r})
            import report_generate
            from pathlib import Path
            report_generate.OUT_DIR = Path({s9_out_dir!r})
            report_generate.FIG_DIR = Path({s9_out_dir!r}) / "figures"
            report_generate.OUT_MD  = Path({s9_out_dir!r}) / "vaccine_report.md"
            report_generate.main()
        """)
        run_script(script, log9, 9)
        st.rerun()

    if st.session_state.step_status[9] == "done":
        report_path = Path(st.session_state.get("s9_out_dir", str(step_default(9)))) / "vaccine_report.md"
        if report_path.exists():
            st.success(f"Report written to `{report_path}`")
            with st.expander("Preview report", expanded=True):
                st.markdown(report_path.read_text(), unsafe_allow_html=False)
