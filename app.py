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

# Widget keys holding per-run paths — cleared on run switch so they re-render
# with correct step_default() values instead of stale paths from the old run.
PER_RUN_PATH_KEYS = [
    "s1_out_dir",
    "s2_tumor_bam", "s2_out_dir",
    "s3_out_dir_manual", "s3_out_dir",
    "s4_vcf", "s4_hla_file", "s4_out_dir",
    "s5_input", "s5_hla_file", "s5_out_dir",
    "s6_input", "s6_hla_file", "s6_out_dir",
    "s7_input", "s7_out_dir",
    "s9_out_dir",
]

def clear_path_widgets():
    for k in PER_RUN_PATH_KEYS:
        st.session_state.pop(k, None)

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
    icons  = {"pending": "–", "running": "···", "done": "✓", "failed": "✕"}
    labels = {"pending": "Pending", "running": "Running", "done": "Complete", "failed": "Failed"}
    st.markdown(f"""
    <div class="step-banner">
        <span class="step-badge">STEP {step_num}</span>
        <span class="step-title">{STEP_NAMES[step_num]}</span>
        <span class="step-status-{status}">{icons[status]} {labels[status]}</span>
    </div>
    """, unsafe_allow_html=True)


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
        clear_path_widgets()
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
                clear_path_widgets()
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
    icons = {"pending": "·", "running": "›", "done": "✓", "failed": "✕"}
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

    st.markdown("""
    <style>
    .pipeline-header { font-size: 1.6rem; font-weight: 700; margin-bottom: 0.2rem; color: #0f1117; }
    .pipeline-sub    { font-size: 0.95rem; color: #555; margin-bottom: 1.5rem; }
    .run-banner {
        background: #f0f4ff; border-left: 4px solid #4a6cf7;
        border-radius: 0 8px 8px 0; padding: 0.7rem 1.1rem;
        margin-bottom: 1.5rem; font-size: 0.85rem; color: #333;
    }
    .run-banner code { background: #dde4ff; padding: 2px 6px; border-radius: 4px; font-size: 0.82rem; }
    .pipe-wrap { display: flex; flex-direction: column; gap: 0; margin: 1rem 0 2rem 0; max-width: 680px; }
    .pipe-input {
        background: #1a1a2e; color: #e0e0ff; border-radius: 10px 10px 0 0;
        padding: 0.7rem 1.2rem; font-size: 0.82rem; font-weight: 600;
        letter-spacing: 0.04em; text-align: center;
    }
    .pipe-arrow { width: 2px; height: 14px; background: #d1d5db; margin: 0 auto; }
    .pipe-step {
        display: flex; align-items: stretch; border: 1px solid #e5e7eb;
        border-radius: 8px; overflow: hidden; margin: 2px 0;
    }
    .pipe-step-num {
        background: #f9fafb; color: #9ca3af; font-size: 0.68rem; font-weight: 700;
        padding: 0.5rem 0.75rem; display: flex; align-items: center;
        justify-content: center; min-width: 52px; border-right: 1px solid #e5e7eb;
        letter-spacing: 0.06em;
    }
    .pipe-step-body { flex: 1; padding: 0.5rem 1rem; display: flex; align-items: center; justify-content: space-between; }
    .pipe-step-name { font-size: 0.88rem; font-weight: 600; color: #111827; }
    .pipe-step-out  { font-size: 0.74rem; color: #6b7280; margin-top: 1px; }
    .pipe-step-tool {
        font-size: 0.68rem; background: #f3f4f6; color: #374151;
        padding: 2px 9px; border-radius: 20px; white-space: nowrap; margin-left: 0.6rem;
    }
    .pipe-output {
        background: #064e3b; color: #a7f3d0; border-radius: 0 0 10px 10px;
        padding: 0.7rem 1.2rem; font-size: 0.82rem; font-weight: 600; text-align: center;
    }
    .status-grid { display: grid; grid-template-columns: repeat(8, 1fr); gap: 8px; margin-top: 0.5rem; }
    .status-card { border: 1px solid #e5e7eb; border-radius: 10px; padding: 0.65rem 0.3rem; text-align: center; background: #fafafa; }
    .status-card.done    { background: #f0fdf4; border-color: #86efac; }
    .status-card.running { background: #fffbeb; border-color: #fcd34d; }
    .status-card.failed  { background: #fff1f2; border-color: #fca5a5; }
    .status-num  { font-size: 0.62rem; font-weight: 700; color: #9ca3af; letter-spacing: 0.08em; text-transform: uppercase; }
    .status-icon { font-size: 0.8rem; font-weight: 700; line-height: 2; letter-spacing: 0.02em; }
    .status-card.done    .status-icon { color: #16a34a; }
    .status-card.running .status-icon { color: #d97706; }
    .status-card.failed  .status-icon { color: #dc2626; }
    .status-card.pending .status-icon { color: #9ca3af; }
    .status-name { font-size: 0.62rem; color: #6b7280; margin-top: 2px; line-height: 1.3; }
    /* Step tab banner */
    .step-banner {
        display: flex; align-items: center; gap: 12px;
        padding: 0.8rem 1rem; background: #f9fafb;
        border: 1px solid #e5e7eb; border-radius: 10px; margin-bottom: 1.2rem;
    }
    .step-badge {
        background: #1a1a2e; color: #e0e0ff; font-size: 0.7rem; font-weight: 700;
        padding: 4px 10px; border-radius: 20px; letter-spacing: 0.06em; white-space: nowrap;
    }
    .step-title { font-size: 1.05rem; font-weight: 600; color: #111827; }
    .step-status-done    { color: #16a34a; font-size: 0.8rem; font-weight: 600; margin-left: auto; }
    .step-status-running { color: #d97706; font-size: 0.8rem; font-weight: 600; margin-left: auto; }
    .step-status-failed  { color: #dc2626; font-size: 0.8rem; font-weight: 600; margin-left: auto; }
    .step-status-pending { color: #9ca3af; font-size: 0.8rem; font-weight: 600; margin-left: auto; }
    </style>
    """, unsafe_allow_html=True)

    st.markdown('<p class="pipeline-header">🧬 mRNA Vaccine Pipeline</p>', unsafe_allow_html=True)
    st.markdown('<p class="pipeline-sub">Tumor sequencing → neoantigen identification → personalized mRNA vaccine construct</p>', unsafe_allow_html=True)

    run_id  = st.session_state.run_id
    run_dir = str(active_run_dir())
    st.markdown(f'<div class="run-banner"><strong>Active run:</strong> <code>{run_id}</code> &nbsp;·&nbsp; <strong>Output:</strong> <code>{run_dir}</code></div>', unsafe_allow_html=True)

    steps_meta = [
        ("1", "Preprocessing",        "Clean BAMs",             "GATK MarkDuplicates"),
        ("2", "Variant Calling",       "Somatic VCF",            "Mutect2 + Filter"),
        ("3", "HLA Typing",            "HLA alleles",            "OptiType"),
        ("4", "Neoantigen Prediction", "Candidate peptides",     "MHCflurry 2.0"),
        ("5", "Immunogenicity Ranking","Ranked candidates",      "IMPROVE weights"),
        ("6", "Epitope Ordering",      "Ordered epitope FASTA",  "Held-Karp TSP"),
        ("7", "mRNA Design",           "Full mRNA construct",    "VaxPress + RNAfold"),
        ("9", "Report",                "vaccine_report.md",      "matplotlib"),
    ]

    html = '<div class="pipe-wrap">'
    html += '<div class="pipe-input">INPUT — Tumor WES + Matched Normal WES</div>'
    for num, name, output, tool in steps_meta:
        html += f'<div class="pipe-arrow"></div>'
        html += f'''<div class="pipe-step">
            <div class="pipe-step-num">STEP {num}</div>
            <div class="pipe-step-body">
                <div><div class="pipe-step-name">{name}</div><div class="pipe-step-out">→ {output}</div></div>
                <span class="pipe-step-tool">{tool}</span>
            </div></div>'''
    html += '<div class="pipe-arrow"></div>'
    html += '<div class="pipe-output">OUTPUT — results/run_&lt;sample&gt;_&lt;timestamp&gt;/</div>'
    html += '</div>'
    st.markdown(html, unsafe_allow_html=True)

    st.markdown("#### Run Status")
    status_map  = st.session_state.step_status
    step_labels = {1:"Preprocess", 2:"Variants", 3:"HLA", 4:"Neoantigens", 5:"Ranking", 6:"Ordering", 7:"mRNA", 9:"Report"}
    icon_map    = {"pending": "—", "running": "···", "done": "✓", "failed": "✕"}

    cards = '<div class="status-grid">'
    for i in ALL_STEPS:
        s = status_map[i]
        cards += f'<div class="status-card {s}"><div class="status-num">Step {i}</div><div class="status-icon">{icon_map[s]}</div><div class="status-name">{step_labels[i]}</div></div>'
    cards += '</div>'
    st.markdown(cards, unsafe_allow_html=True)

# ─── Step 1: Preprocessing ────────────────────────────────────────────────────

if not entry_b:
    with tabs[tab_idx]:
        tab_idx += 1
        step_header(1)
        st.caption("Convert raw sequencing reads into clean, analysis-ready BAMs using GATK MarkDuplicates and SortSam.")

        # Primary — always visible
        s1_input_bam = st.text_input(
            "Input BAM",
            str(DATA_DIR / "test" / "tumor_chr17.bam"),
            key="s1_input_bam",
        )

        # Advanced — collapsed
        with st.expander("Advanced", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                s1_out_dir = st.text_input("Output directory", str(step_default(1)), key="s1_out_dir")
            with col2:
                s1_gatk_jar  = st.text_input("GATK JAR", gatk_jar, key="s1_gatk_jar")
                s1_reference = st.text_input("Reference (b37)", ref_b37, key="s1_reference")

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
        st.caption("Identify somatic mutations present in tumor but not normal using GATK Mutect2 and FilterMutectCalls.")

        # Primary — always visible
        col1, col2 = st.columns(2)
        with col1:
            s2_tumor_bam = st.text_input("Tumor BAM", str(step_default(1) / "sorted.bam"), key="s2_tumor_bam")
        with col2:
            s2_normal_bam = st.text_input("Normal BAM", str(DATA_DIR / "test" / "normal_chr17.bam"), key="s2_normal_bam")

        # Parameters — collapsed
        with st.expander("Parameters", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                s2_tumor_name  = st.text_input("Tumor sample name",  "HCC1143_tumor",  key="s2_tumor_name")
                s2_normal_name = st.text_input("Normal sample name", "HCC1143_normal", key="s2_normal_name")
            with col2:
                s2_intervals = st.text_input("Intervals (chromosome)", "chr17", key="s2_intervals")

        # Advanced — collapsed
        with st.expander("Advanced", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                s2_out_dir   = st.text_input("Output directory", str(step_default(2)), key="s2_out_dir")
            with col2:
                s2_gatk_jar  = st.text_input("GATK JAR", gatk_jar, key="s2_gatk_jar")
                s2_reference = st.text_input("Reference (hg38)", ref_hg38, key="s2_reference")

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
    st.caption("Enter HLA alleles directly if you have them from clinical typing, or run OptiType via CLI.")

    # Primary — manual entry is the main path
    st.markdown("**HLA alleles**")
    manual_hla = st.text_area(
        "One allele per line (e.g. HLA-A*02:01)",
        "HLA-A*02:01\nHLA-A*24:02\nHLA-B*07:02\nHLA-B*35:01\nHLA-C*04:01\nHLA-C*07:02",
        key="manual_hla",
        height=140,
    )

    # Advanced — output dir
    with st.expander("Advanced", expanded=False):
        s3_out_dir_manual = st.text_input(
            "Output directory", str(step_default(3)), key="s3_out_dir_manual"
        )

    if st.button("Write HLA alleles", type="primary"):
        out = Path(s3_out_dir_manual)
        out.mkdir(parents=True, exist_ok=True)
        alleles = [a.strip() for a in manual_hla.strip().splitlines() if a.strip()]
        with open(out / "hla_alleles.txt", "w") as f:
            f.write("\n".join(alleles) + "\n")
        st.success(f"Wrote {len(alleles)} alleles to {out / 'hla_alleles.txt'}")
        st.session_state.step_status[3] = "done"

    # OptiType fallback — fully collapsed
    with st.expander("OptiType fallback (CLI — no HLA from clinical typing)", expanded=False):
        st.info(
            "```bash\n"
            "# Phase 1 — extract HLA reads\n"
            "docker compose run --rm pipeline python3 scripts/hla_typing.py --optitype\n\n"
            "# Phase 2 — run OptiType container\n"
            "OPTITYPE_DATA_DIR=./results/run_<id>/step3 \\\n"
            "OPTITYPE_SAMPLE=hcc1143_normal \\\n"
            "docker compose run --rm optitype\n\n"
            "# Phase 3 — parse results\n"
            "docker compose run --rm pipeline python3 scripts/hla_typing.py --parse\n"
            "```"
        )
        col1, col2 = st.columns(2)
        with col1:
            s3_out_dir = st.text_input("Step 3 output directory", str(step_default(3)), key="s3_out_dir")
        with col2:
            s3_sample_prefix = st.text_input("Sample prefix", "hcc1143_normal", key="s3_sample_prefix")

        show_prior_log(3)
        log3 = st.empty()
        if st.button("Parse OptiType results → hla_alleles.txt"):
            script = textwrap.dedent(f"""\
                import sys
                sys.path.insert(0, {str(SCRIPTS_DIR)!r})
                import hla_typing
                from pathlib import Path
                out_dir = Path({s3_out_dir!r})
                out_dir.mkdir(parents=True, exist_ok=True)
                print("Step 3 — parsing OptiType results")
                hla_alleles = hla_typing.parse_optitype(out_dir, {s3_sample_prefix!r})
                alleles_file = out_dir / "hla_alleles.txt"
                with open(alleles_file, "w") as f:
                    for allele in hla_alleles:
                        f.write(allele + "\\n")
                print("── HLA Alleles ──")
                for allele in hla_alleles:
                    print(f"  {{allele}}")
                print(f"Step 3 complete. Alleles saved to: {{alleles_file}}")
            """)
            run_script(script, log3, 3)
            st.rerun()

# ─── Step 4: Neoantigen Prediction ────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(4)
    st.caption("Extract mutant peptides from somatic variants and predict MHC binding affinity with MHCflurry 2.0.")

    # Primary — inputs
    col1, col2 = st.columns(2)
    with col1:
        s4_vcf = st.text_input(
            "VCF / MAF file",
            str(step_default(2) / "filtered_variants.vcf.gz"),
            key="s4_vcf",
        )
    with col2:
        s4_hla_file = st.text_input(
            "HLA alleles file",
            str(step_default(3) / "hla_alleles.txt"),
            key="s4_hla_file",
        )

    # Parameters — collapsed
    with st.expander("Parameters", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            s4_peptide_lengths = st.multiselect(
                "Peptide lengths",
                options=[8, 9, 10, 11, 12],
                default=[8, 9, 10, 11],
                key="s4_peptide_lengths",
            )
        with col2:
            s4_binding_cutoff = st.number_input(
                "Binding cutoff (IC50 nM)",
                min_value=50, max_value=5000, value=500, step=50,
                key="s4_binding_cutoff",
                help="MHC binding threshold — standard is 500 nM",
            )

    # Advanced — collapsed
    with st.expander("Advanced", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            s4_out_dir = st.text_input("Output directory", str(step_default(4)), key="s4_out_dir")
        with col2:
            s4_reference = st.text_input("Reference (hg38)", ref_hg38, key="s4_reference")

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
    st.caption("Score candidates using a composite immunogenicity model based on IMPROVE 2024 feature importance.")

    # Primary — input + top N
    col1, col2 = st.columns([3, 1])
    with col1:
        s5_input = st.text_input(
            "Candidates TSV (Step 4 output)",
            str(step_default(4) / "candidate_neoantigens.tsv"),
            key="s5_input",
        )
    with col2:
        s5_top_n = st.number_input("Top N", min_value=5, max_value=50, value=30, step=1, key="s5_top_n")

    # Parameters — scoring weights, open by default since they're interesting
    with st.expander("Scoring weights", expanded=True):
        st.caption("Based on IMPROVE (Frontiers Immunology 2024) feature importance. Must sum to 1.0.")
        col1, col2 = st.columns(2)
        with col1:
            s5_w_pres = st.slider("Presentation score", 0.0, 1.0, 0.40, 0.05, key="s5_w_pres")
            s5_w_agr  = st.slider("Agretopicity (DAI)", 0.0, 1.0, 0.25, 0.05, key="s5_w_agr")
            s5_w_vaf  = st.slider("VAF / clonality",    0.0, 1.0, 0.20, 0.05, key="s5_w_vaf")
        with col2:
            s5_w_bl   = st.slider("BLOSUM score",       0.0, 1.0, 0.10, 0.05, key="s5_w_bl")
            s5_w_for  = st.slider("Foreignness",        0.0, 1.0, 0.05, 0.05, key="s5_w_for")
            total_w = s5_w_pres + s5_w_agr + s5_w_vaf + s5_w_bl + s5_w_for
            color = "green" if abs(total_w - 1.0) < 0.01 else "red"
            st.markdown(f"Total: <span style='color:{color}; font-weight:600'>{total_w:.2f}</span>", unsafe_allow_html=True)

    # Advanced — filters, output dir, HLA file
    with st.expander("Advanced", expanded=False):
        col1, col2, col3 = st.columns(3)
        with col1:
            s5_out_dir = st.text_input("Output directory", str(step_default(5)), key="s5_out_dir")
        with col2:
            s5_hla_file = st.text_input("HLA alleles file", str(step_default(3) / "hla_alleles.txt"), key="s5_hla_file")
        with col3:
            s5_max_aff  = st.number_input("Max affinity (nM)", 50, 5000, 500, 50, key="s5_max_aff")
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
    st.caption("Optimal epitope concatenation order to minimize junctional neoantigens. Held-Karp exact TSP for N≤15, greedy for N>15.")

    # Primary — input file
    s6_input = st.text_input(
        "Ranked neoantigens TSV (Step 5 output)",
        str(step_default(5) / "ranked_neoantigens.tsv"),
        key="s6_input",
    )

    # Parameters — linker config
    with st.expander("Parameters", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            s6_linker       = st.text_input("GPGPG linker sequence", "GPGPG", key="s6_linker")
        with col2:
            s6_junction_len = st.number_input(
                "Junction peptide length", min_value=7, max_value=11, value=9, step=1,
                key="s6_junction_len", help="Length of junctional peptides to check for MHC binding",
            )
            s6_max_greedy = st.number_input(
                "Exact TSP threshold (N)",  min_value=5, max_value=50, value=15, step=1,
                key="s6_max_greedy", help="Use exact Held-Karp below this, greedy above",
            )

    # Advanced — paths
    with st.expander("Advanced", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            s6_out_dir  = st.text_input("Output directory", str(step_default(6)), key="s6_out_dir")
        with col2:
            s6_hla_file = st.text_input("HLA alleles file", str(step_default(3) / "hla_alleles.txt"), key="s6_hla_file")

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
            epitope_ordering.MAX_EXACT     = {s6_max_greedy}
            epitope_ordering.main()
        """)
        run_script(script, log6, 6)
        st.rerun()

# ─── Step 7: mRNA Design ──────────────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(7)
    st.caption("VaxPress codon optimization across multiple profiles, assembled with 5'UTR/Kozak/poly-A and validated with RNAfold.")

    # Primary — input FASTA
    s7_input = st.text_input(
        "Ordered epitopes FASTA (Step 6 output)",
        str(step_default(6) / "ordered_epitopes.fasta"),
        key="s7_input",
    )

    # Candidate profiles — open by default, this is the interesting part
    with st.expander("Candidate profiles", expanded=True):
        st.caption("Each profile sweeps the MFE/CAI tradeoff (mirroring LinearDesign Nature 2023 λ sweep).")
        s7_candidates = []
        profile_defaults = [
            ("A_balanced",    "Balanced MFE + CAI", {}),
            ("B_stability",   "Stability-biased — longer half-life", {"mfe": 10, "cai": 2, "loop": 5}),
            ("C_expression",  "Expression-biased — higher protein yield", {"mfe": 2, "cai": 10, "gc": 3}),
            ("D_low_uridine", "Uridine-minimized — reduced innate immune activation", {"ucount": 10, "mfe": 5, "cai": 3}),
        ]
        for i, (pid, pdesc, pweights) in enumerate(profile_defaults):
            with st.container():
                cols = st.columns([1, 3, 4])
                with cols[0]:
                    enabled = st.checkbox("On", value=True, key=f"s7_cand_enabled_{i}")
                    cand_id = st.text_input("ID", pid, key=f"s7_cand_id_{i}", label_visibility="collapsed")
                with cols[1]:
                    cand_desc = st.text_input("Description", pdesc, key=f"s7_cand_desc_{i}", label_visibility="collapsed")
                with cols[2]:
                    weight_keys = ["mfe", "cai", "ucount", "loop", "gc"]
                    weight_vals = {}
                    wcols = st.columns(5)
                    for j, wk in enumerate(weight_keys):
                        with wcols[j]:
                            default_val = pweights.get(wk, 0)
                            v = st.number_input(wk, min_value=0, max_value=20, value=default_val, step=1, key=f"s7_w_{i}_{wk}")
                            if v > 0:
                                weight_vals[wk] = v
                if enabled:
                    s7_candidates.append({"id": cand_id, "description": cand_desc, "weights": weight_vals})
                st.divider()

    # Parameters — VaxPress tuning
    with st.expander("Parameters", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            s7_vaxpress_iters = st.number_input("VaxPress iterations", min_value=50, max_value=1000, value=200, step=50, key="s7_vaxpress_iters")
        with col2:
            s7_vaxpress_procs = st.number_input("VaxPress processes",  min_value=1,  max_value=32,   value=8,   step=1,  key="s7_vaxpress_procs")

    # Advanced — construct sequences, output dir
    with st.expander("Advanced", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            s7_out_dir    = st.text_input("Output directory", str(step_default(7)), key="s7_out_dir")
            s7_poly_a_len = st.number_input("Poly-A tail length (nt)", min_value=60, max_value=300, value=120, step=10, key="s7_poly_a_len")
        with col2:
            s7_utr5  = st.text_input("5' UTR",         "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC", key="s7_utr5")
            s7_kozak = st.text_input("Kozak sequence",  "GCCACCATG", key="s7_kozak")

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

    if st.session_state.step_status[7] == "done":
        st.subheader("Results")
        out_dir = Path(st.session_state.get("s7_out_dir", str(step_default(7))))
        comparison_file = out_dir / "candidate_comparison.json"
        if comparison_file.exists():
            with open(comparison_file) as f:
                comparison = json.load(f)
            import pandas as pd
            rows = []
            for m in comparison:
                rows.append({
                    "Candidate":      m["candidate_id"],
                    "Description":    m["description"],
                    "GC %":           f"{m['gc_content']:.1%}" if m.get("gc_content") else "N/A",
                    "MFE (kcal/mol)": f"{m['mfe_kcal_mol']:.2f}" if m.get("mfe_kcal_mol") else "N/A",
                    "Length (nt)":    m.get("total_length_nt", "N/A"),
                    "QC":             "ISSUES" if m.get("issues") else ("warn" if m.get("warnings") else "pass"),
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True)

# ─── Step 9: Report Generation ────────────────────────────────────────────────

with tabs[tab_idx]:
    tab_idx += 1
    step_header(9)
    st.caption("Generate a Markdown report with dynamic tables and matplotlib plots from all pipeline outputs.")

    with st.expander("Advanced", expanded=False):
        s9_out_dir = st.text_input("Output directory", str(step_default(9)), key="s9_out_dir")

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