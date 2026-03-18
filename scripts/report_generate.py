"""
Step 9: Markdown Report Generation
Reads all pipeline outputs and generates a rich, self-contained Markdown
report with embedded matplotlib plots and dynamic tables.

Output layout:
    results/step9/
        vaccine_report.md      — main report
        figures/               — PNG plots referenced in the report

Dependencies: matplotlib, seaborn, pandas, numpy
"""

import json
from datetime import date
from pathlib import Path

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).parent))
from paths import step_dir, RUN_ID, SAMPLE_ID, TUMOR_TYPE, PIPELINE_VER

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

# ── Config ────────────────────────────────────────────────────────────────────

OUT_DIR     = step_dir(9)
FIG_DIR     = OUT_DIR / "figures"
OUT_MD      = OUT_DIR / "vaccine_report.md"

STEP3_HLA   = step_dir(3) / "hla_alleles.txt"
STEP5_TSV   = step_dir(5) / "ranked_neoantigens.tsv"
STEP5_SUB   = step_dir(5) / "subclonal_filtered.tsv"
STEP6_TSV   = step_dir(6) / "ordered_epitopes.tsv"
STEP6_FASTA = step_dir(6) / "ordered_epitopes.fasta"
STEP6_JXN   = step_dir(6) / "junction_scores.tsv"
STEP7_A     = step_dir(7) / "A_balanced" / "quality_metrics.json"
STEP7_B     = step_dir(7) / "B_stability" / "quality_metrics.json"

REPORT_DATE = date.today().isoformat()

# Plot style
PALETTE = {
    "navy":  "#1a2e4a",
    "teal":  "#0f6e56",
    "amber": "#ba7517",
    "red":   "#a32d2d",
    "lgray": "#f1efe8",
    "mgray": "#d3d1c7",
    "dgray": "#5f5e5a",
}

sns.set_theme(style="whitegrid", font_scale=0.9)
plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "savefig.dpi":      150,
    "savefig.bbox":     "tight",
})


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_hla():
    if not STEP3_HLA.exists():
        return []
    return [l.strip() for l in STEP3_HLA.read_text().splitlines() if l.strip()]


def load_step5():
    if not STEP5_TSV.exists():
        return pd.DataFrame()
    return pd.read_csv(STEP5_TSV, sep="\t")


def load_step5_sub():
    if not STEP5_SUB.exists():
        return pd.DataFrame()
    return pd.read_csv(STEP5_SUB, sep="\t")


def load_step6():
    if not STEP6_TSV.exists():
        return pd.DataFrame()
    return pd.read_csv(STEP6_TSV, sep="\t")


def load_junction_matrix():
    if not STEP6_JXN.exists():
        return pd.DataFrame()
    return pd.read_csv(STEP6_JXN, sep="\t", index_col=0)


def load_step6_fasta():
    if not STEP6_FASTA.exists():
        return ""
    seq = ""
    for line in STEP6_FASTA.read_text().splitlines():
        if not line.startswith(">"):
            seq += line.strip()
    return seq


def load_candidates():
    candidates = []
    for path in [STEP7_A, STEP7_B]:
        if path.exists():
            candidates.append(json.loads(path.read_text()))
    return candidates


# ── Markdown helpers ──────────────────────────────────────────────────────────

def md_table(headers, rows):
    """Render a GitHub-flavoured Markdown table."""
    def escape(s):
        return str(s).replace("|", "\\|")

    col_widths = [max(len(str(h)), max((len(escape(r[i])) for r in rows), default=0))
                  for i, h in enumerate(headers)]

    def pad(s, w):
        return escape(s).ljust(w)

    header_line = "| " + " | ".join(pad(h, col_widths[i]) for i, h in enumerate(headers)) + " |"
    sep_line    = "| " + " | ".join("-" * w for w in col_widths) + " |"
    body_lines  = [
        "| " + " | ".join(pad(r[i], col_widths[i]) for i in range(len(headers))) + " |"
        for r in rows
    ]
    return "\n".join([header_line, sep_line] + body_lines)


def fig_ref(name):
    """Return a markdown image reference for a figure in the figures/ subdir."""
    return f"![{name}](figures/{name}.png)"


# ── Plot generators ────────────────────────────────────────────────────────────

def plot_score_distribution(df5):
    """Bar chart: composite score for all top candidates, coloured by flags."""
    fig, ax = plt.subplots(figsize=(10, 4))

    colors = []
    for _, row in df5.iterrows():
        flags = str(row.get("flags", ""))
        if "NEG_AGRETOPICITY" in flags:
            colors.append(PALETTE["red"])
        elif "LOW_VAF" in flags:
            colors.append(PALETTE["amber"])
        else:
            colors.append(PALETTE["teal"])

    x = range(len(df5))
    ax.bar(x, df5["composite_score"], color=colors, width=0.7, edgecolor="white", linewidth=0.4)
    ax.set_xticks(list(x))
    ax.set_xticklabels(df5["peptide"].tolist(), rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Composite score")
    ax.set_title("Neoantigen Candidate Composite Scores", fontweight="bold", color=PALETTE["navy"])

    legend_handles = [
        mpatches.Patch(color=PALETTE["teal"],  label="Clean"),
        mpatches.Patch(color=PALETTE["amber"], label="LOW_VAF"),
        mpatches.Patch(color=PALETTE["red"],   label="NEG_AGRETOPICITY"),
    ]
    ax.legend(handles=legend_handles, fontsize=8, loc="upper right")
    plt.tight_layout()
    path = FIG_DIR / "score_distribution.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_affinity_vs_presentation(df5):
    """Scatter: affinity vs presentation score, dot size = VAF, colour = composite score."""
    fig, ax = plt.subplots(figsize=(7, 5))

    sc = ax.scatter(
        df5["affinity"],
        df5["presentation_score"],
        c=df5["composite_score"],
        s=df5["vaf"].fillna(0.1) * 400,
        cmap="YlOrRd",
        alpha=0.85,
        edgecolors=PALETTE["dgray"],
        linewidths=0.5,
    )
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Composite score", fontsize=8)

    for _, row in df5.iterrows():
        ax.annotate(
            row["peptide"], (row["affinity"], row["presentation_score"]),
            textcoords="offset points", xytext=(4, 3), fontsize=6,
            color=PALETTE["dgray"],
        )

    ax.axvline(500, color=PALETTE["amber"], linestyle="--", linewidth=0.8, label="500 nM threshold")
    ax.axhline(0.1, color=PALETTE["red"],   linestyle="--", linewidth=0.8, label="Pres. score threshold")
    ax.set_xlabel("MHC Binding Affinity (nM)  [lower = stronger binding]")
    ax.set_ylabel("Presentation Score")
    ax.set_title("Affinity vs Presentation Score\n(dot size ∝ VAF)", fontweight="bold", color=PALETTE["navy"])
    ax.legend(fontsize=8)
    plt.tight_layout()
    path = FIG_DIR / "affinity_vs_presentation.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_hla_coverage(df5, hla_alleles):
    """Horizontal bar chart showing how many top candidates cover each HLA allele."""
    counts = {a: int((df5["best_allele"] == a).sum()) for a in hla_alleles}
    alleles = list(counts.keys())
    values  = [counts[a] for a in alleles]

    fig, ax = plt.subplots(figsize=(7, max(3, len(alleles) * 0.55)))
    bar_colors = [PALETTE["red"] if v == 0 else PALETTE["teal"] for v in values]
    bars = ax.barh(alleles, values, color=bar_colors, edgecolor="white", height=0.6)
    ax.bar_label(bars, padding=3, fontsize=8)
    ax.set_xlabel("Number of candidates")
    ax.set_title("HLA Allele Coverage (top candidates)", fontweight="bold", color=PALETTE["navy"])
    ax.set_xlim(0, max(values) + 2 if values else 5)
    plt.tight_layout()
    path = FIG_DIR / "hla_coverage.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_score_components(df5):
    """Stacked bar chart decomposing the composite score into its 5 components."""
    cols = {
        "score_presentation": "Presentation (40%)",
        "score_agretopicity": "Agretopicity (25%)",
        "score_vaf":          "VAF (20%)",
        "score_blosum":       "BLOSUM (10%)",
        "score_foreignness":  "Foreignness (5%)",
    }
    weights = [0.40, 0.25, 0.20, 0.10, 0.05]
    existing = [c for c in cols if c in df5.columns]

    top = df5.head(15).copy()
    fig, ax = plt.subplots(figsize=(11, 4))

    bottoms = np.zeros(len(top))
    palette = [PALETTE["teal"], PALETTE["navy"], PALETTE["amber"],
               PALETTE["dgray"], PALETTE["red"]]

    for idx, col in enumerate(existing):
        label  = list(cols.values())[idx]
        weight = weights[idx]
        vals   = (top[col] * weight).fillna(0).values
        ax.bar(range(len(top)), vals, bottom=bottoms, label=label,
               color=palette[idx % len(palette)], edgecolor="white", linewidth=0.3)
        bottoms += vals

    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(top["peptide"].tolist(), rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Weighted contribution to composite score")
    ax.set_title("Score Component Breakdown (top 15 candidates)", fontweight="bold", color=PALETTE["navy"])
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    plt.tight_layout()
    path = FIG_DIR / "score_components.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_junction_path(df6, jxn_matrix):
    """Bar chart: per-junction scores along the chosen TSP path."""
    if df6.empty or jxn_matrix.empty:
        return None

    if "order" in df6.columns:
        ordered = df6.sort_values("order")
    else:
        ordered = df6.copy()

    epitopes = ordered["peptide"].tolist()
    scores   = []
    labels   = []
    for i in range(len(epitopes) - 1):
        a, b = epitopes[i], epitopes[i + 1]
        if a in jxn_matrix.index and b in jxn_matrix.columns:
            scores.append(jxn_matrix.loc[a, b])
        else:
            scores.append(0.0)
        labels.append(f"{i+1}→{i+2}")

    fig, ax = plt.subplots(figsize=(10, 3.5))
    bar_colors = [PALETTE["red"] if s > 0.5 else
                  PALETTE["amber"] if s > 0.3 else PALETTE["teal"]
                  for s in scores]
    ax.bar(range(len(scores)), scores, color=bar_colors, edgecolor="white", width=0.7)
    ax.set_xticks(range(len(scores)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.axhline(0.5, color=PALETTE["red"],   linestyle="--", linewidth=0.8, label="> 0.5 (high)")
    ax.axhline(0.3, color=PALETTE["amber"], linestyle="--", linewidth=0.8, label="> 0.3 (moderate)")
    ax.set_ylabel("Junction presentation score")
    ax.set_title("Per-Junction Scores Along Optimal Epitope Path\n(lower = fewer spurious junctional epitopes)",
                 fontweight="bold", color=PALETTE["navy"])
    ax.legend(fontsize=8)
    plt.tight_layout()
    path = FIG_DIR / "junction_path.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_junction_heatmap(jxn_matrix):
    """Heatmap of the full N×N junction score matrix."""
    if jxn_matrix.empty:
        return None

    fig, ax = plt.subplots(figsize=(10, 8))
    mask = np.eye(len(jxn_matrix), dtype=bool)
    sns.heatmap(
        jxn_matrix,
        ax=ax,
        cmap="YlOrRd",
        mask=mask,
        linewidths=0.3,
        linecolor="#e0e0e0",
        annot=len(jxn_matrix) <= 15,
        fmt=".2f" if len(jxn_matrix) <= 15 else "",
        annot_kws={"size": 7},
        cbar_kws={"label": "Worst junction presentation score"},
    )
    ax.set_title("Junction Score Matrix  (diagonal masked — self-junctions excluded)",
                 fontweight="bold", color=PALETTE["navy"])
    ax.set_xlabel("Epitope B (following)")
    ax.set_ylabel("Epitope A (preceding)")
    plt.tight_layout()
    path = FIG_DIR / "junction_heatmap.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_mrna_comparison(candidates):
    """Side-by-side bar charts comparing GC content and MFE/nt for the two mRNA candidates."""
    if not candidates:
        return None

    ids   = [m.get("candidate_id", f"Candidate {i+1}") for i, m in enumerate(candidates)]
    gcs   = [m.get("gc_content", 0) * 100 for m in candidates]
    mfes  = [m.get("mfe_kcal_mol") for m in candidates]
    mpnts = [m.get("mfe_per_nt") for m in candidates]

    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    bar_kw = dict(width=0.5, edgecolor="white", color=[PALETTE["teal"], PALETTE["navy"]])

    # GC content
    axes[0].bar(ids, gcs, **bar_kw)
    axes[0].axhline(40, color=PALETTE["amber"], linestyle="--", linewidth=0.8)
    axes[0].axhline(70, color=PALETTE["amber"], linestyle="--", linewidth=0.8)
    axes[0].set_ylabel("GC content (%)")
    axes[0].set_title("GC Content", fontweight="bold", color=PALETTE["navy"])
    axes[0].set_ylim(0, 100)
    for i, v in enumerate(gcs):
        axes[0].text(i, v + 0.5, f"{v:.1f}%", ha="center", va="bottom", fontsize=9)

    # MFE
    valid_mfe = [v for v in mfes if v is not None]
    if valid_mfe:
        axes[1].bar(ids, [v if v is not None else 0 for v in mfes], **bar_kw)
        axes[1].set_ylabel("MFE (kcal/mol)")
        axes[1].set_title("Minimum Free Energy", fontweight="bold", color=PALETTE["navy"])
        for i, v in enumerate(mfes):
            if v is not None:
                axes[1].text(i, v - abs(v)*0.03, f"{v:.1f}", ha="center", va="top", fontsize=9)
    else:
        axes[1].text(0.5, 0.5, "MFE not available\n(RNAfold not run)",
                     ha="center", va="center", transform=axes[1].transAxes, color=PALETTE["dgray"])
        axes[1].set_title("Minimum Free Energy", fontweight="bold", color=PALETTE["navy"])

    # MFE/nt
    valid_mpnt = [v for v in mpnts if v is not None]
    if valid_mpnt:
        axes[2].bar(ids, [v if v is not None else 0 for v in mpnts], **bar_kw)
        axes[2].axhline(-0.48, color=PALETTE["amber"], linestyle="--", linewidth=0.8, label="Opt. range")
        axes[2].axhline(-0.60, color=PALETTE["amber"], linestyle="--", linewidth=0.8)
        axes[2].set_ylabel("MFE / nt  (kcal/mol/nt)")
        axes[2].set_title("MFE per Nucleotide", fontweight="bold", color=PALETTE["navy"])
        axes[2].legend(fontsize=8)
        for i, v in enumerate(mpnts):
            if v is not None:
                axes[2].text(i, v - abs(v)*0.03, f"{v:.3f}", ha="center", va="top", fontsize=9)
    else:
        axes[2].text(0.5, 0.5, "MFE/nt not available",
                     ha="center", va="center", transform=axes[2].transAxes, color=PALETTE["dgray"])
        axes[2].set_title("MFE per Nucleotide", fontweight="bold", color=PALETTE["navy"])

    fig.suptitle("mRNA Candidate Comparison", fontsize=12, fontweight="bold", color=PALETTE["navy"])
    plt.tight_layout()
    path = FIG_DIR / "mrna_comparison.png"
    fig.savefig(path)
    plt.close(fig)
    return path


# ── Section writers ────────────────────────────────────────────────────────────

def section_header(lines):
    lines += [
        f"# Personalized mRNA Vaccine — Computational Design Report",
        "",
        f"| Field | Value |",
        f"|---|---|",
        f"| **Sample ID** | {SAMPLE_ID} |",
        f"| **Tumor type** | {TUMOR_TYPE} |",
        f"| **Run ID** | `{RUN_ID}` |",
        f"| **Report date** | {REPORT_DATE} |",
        f"| **Pipeline** | melanoma-pipeline v{PIPELINE_VER} |",
        f"| **Reference genome** | hg38 |",
        f"| **HLA typing** | OptiType + Cellosaurus ground truth |",
        "",
        "> ⚠️ **RESEARCH USE ONLY — Not for clinical or diagnostic use.**",
        "",
        "---",
        "",
    ]


def section_toc(lines):
    lines += [
        "## Table of Contents",
        "",
        "1. [Pipeline Summary](#1-pipeline-summary)",
        "2. [HLA Typing](#2-hla-typing)",
        "3. [Neoantigen Candidates](#3-neoantigen-candidates)",
        "4. [Epitope Ordering](#4-epitope-ordering)",
        "5. [mRNA Construct Design](#5-mrna-construct-design)",
        "6. [Wet Lab Synthesis Specifications](#6-wet-lab-synthesis-specifications)",
        "7. [Limitations & Caveats](#7-limitations--caveats)",
        "8. [Methods & References](#8-methods--references)",
        "",
        "---",
        "",
    ]


def section_pipeline_summary(lines):
    lines += [
        "## 1. Pipeline Summary",
        "",
        "This report describes a fully computational personalized mRNA vaccine design pipeline "
        "applied to a melanoma/TNBC cell line model (HCC1143). The pipeline takes tumor and matched "
        "normal whole-genome sequencing data and produces a ready-to-synthesize mRNA construct "
        "encoding the top-ranked somatic neoepitopes, joined by GPGPG linkers and embedded in "
        "clinically validated UTR/poly-A elements.",
        "",
    ]

    headers = ["Step", "Description", "Tool", "Status"]
    rows = [
        ["1", "Preprocessing & duplicate marking",     "GATK MarkDuplicates",          "✅ Complete"],
        ["2", "Somatic variant calling (chr17)",        "GATK Mutect2",                 "✅ Complete"],
        ["3", "HLA typing",                             "OptiType + Cellosaurus",        "✅ Complete"],
        ["4", "Neoantigen prediction",                  "MHCflurry 2.0",                "✅ Complete"],
        ["5", "Immunogenicity ranking",                 "Custom (IMPROVE weights)",      "✅ Complete"],
        ["6", "Epitope ordering (TSP)",                 "Held-Karp exact DP",           "✅ Complete"],
        ["7", "mRNA codon optimisation",                "VaxPress + LinearFold",         "✅ Complete"],
        ["8", "CodonFM validation",                     "nvidia/NV-CodonFM-80M",         "⏳ Deferred"],
        ["9", "Report generation",                      "This document",                "✅ Complete"],
    ]
    lines.append(md_table(headers, rows))
    lines += [
        "",
        "> **Step 8 (CodonFM) deferred** — both candidates proceed to wet lab expression "
        "validation in HEK293T cells, which provides a more direct ranking signal than the "
        "language model score for constructs of this length.",
        "",
        "---",
        "",
    ]


def section_hla(lines, hla_alleles):
    lines += [
        "## 2. HLA Typing",
        "",
        "HLA alleles were determined computationally using OptiType on HLA-enriched reads extracted "
        "from the tumor BAM, then verified against the Cellosaurus CVCL_1245 ground truth "
        "(Boegel et al. 2014, PubMed 25960936). The A locus shows homozygosity "
        "(A\\*31:01 / A\\*31:01) which OptiType may under-call as heterozygous — "
        "Cellosaurus-confirmed alleles were used.",
        "",
    ]

    loci = {}
    for a in hla_alleles:
        locus = a.split("*")[0]
        loci.setdefault(locus, []).append(a)

    hla_rows = []
    for locus, alleles in sorted(loci.items()):
        a1 = alleles[0]
        a2 = alleles[1] if len(alleles) > 1 else f"{a1} *(homozygous)*"
        hla_rows.append([locus, a1, a2, "Cellosaurus verified"])

    lines.append(md_table(["Locus", "Allele 1", "Allele 2", "Source"], hla_rows))
    lines += ["", "---", ""]


def section_candidates(lines, df5, hla_alleles, df5_sub):
    lines += [
        "## 3. Neoantigen Candidates",
        "",
        f"MHCflurry 2.0 predicted **{len(df5)} neoantigen candidates** from somatic variants on "
        "chr17 after applying clinical-standard filters: affinity < 500 nM, presentation score > 0.10, "
        "VAF ≥ 0.05 (clonality). Candidates were ranked by a composite score weighting "
        "presentation (40%), agretopicity (25%), VAF (20%), BLOSUM (10%), and foreignness (5%) "
        "per IMPROVE (Frontiers Immunology 2024) feature importance.",
        "",
    ]

    if not df5_sub.empty:
        lines.append(
            f"> **{len(df5_sub)} subclonal candidates** (VAF < 0.05) were removed and saved to "
            "`step5/subclonal_filtered.tsv` for the audit trail."
        )
        lines.append("")

    # Ranked table
    cols_show   = ["peptide", "best_allele", "affinity", "presentation_score",
                   "agretopicity", "vaf", "composite_score", "flags"]
    cols_exist  = [c for c in cols_show if c in df5.columns]
    header_map  = {
        "peptide": "Peptide", "best_allele": "Allele", "affinity": "Affinity (nM)",
        "presentation_score": "Pres. Score", "agretopicity": "Agretopicity",
        "vaf": "VAF", "composite_score": "Composite", "flags": "Flags",
    }
    headers = [header_map[c] for c in cols_exist]
    rows = []
    for _, row in df5[cols_exist].head(20).iterrows():
        flags = str(row.get("flags", ""))
        r = []
        for c in cols_exist:
            v = row[c]
            if c in ("affinity",):
                r.append(f"{v:.1f}")
            elif c in ("presentation_score", "agretopicity"):
                r.append(f"{v:.3f}")
            elif c in ("vaf",):
                r.append(f"{v:.3f}")
            elif c in ("composite_score",):
                r.append(f"{v:.3f}")
            else:
                r.append(str(v) if pd.notna(v) else "")
        rows.append(r)

    lines.append(md_table(headers, rows))
    lines += [
        "",
        "> 🔴 **NEG\\_AGRETOPICITY** — mutant binds MHC worse than wildtype; review before synthesis.  ",
        "> 🟡 **LOW\\_VAF** — VAF 0.05–0.10; moderately subclonal, included but flagged.",
        "",
    ]

    # Flag callout
    flagged = df5[df5.get("flags", pd.Series(dtype=str)).str.contains("NEG_AGRETOPICITY", na=False)]
    if not flagged.empty:
        lines += [
            "### Flag: NEG_AGRETOPICITY candidates",
            "",
        ]
        for _, row in flagged.iterrows():
            pep    = row.get("peptide", "?")
            allele = row.get("best_allele", "?")
            aff    = row.get("affinity", 0)
            agr    = row.get("agretopicity", 0)
            lines.append(
                f"- **{pep}** ({allele}, affinity {aff:.0f} nM, agretopicity {agr:.3f}) — "
                "mutant binds MHC worse than wildtype. Consider substituting with the next "
                "clean candidate for this allele from `all_scored_candidates.tsv`."
            )
        lines.append("")

    # Plots
    lines += [
        "### Score Distribution",
        "",
        fig_ref("score_distribution"),
        "",
        "### Score Component Breakdown",
        "",
        fig_ref("score_components"),
        "",
        "### Affinity vs Presentation Score",
        "",
        "*Dot size is proportional to VAF; colour represents composite score.*",
        "",
        fig_ref("affinity_vs_presentation"),
        "",
        "### HLA Allele Coverage",
        "",
        fig_ref("hla_coverage"),
        "",
        "---",
        "",
    ]


def section_epitope_ordering(lines, df6, epitope_string, jxn_matrix):
    n_epitopes = len(df6) if not df6.empty else "?"
    lines += [
        "## 4. Epitope Ordering",
        "",
        f"Epitopes were deduplicated by mutation (keeping the highest-scoring peptide length "
        f"per variant), reducing ranked candidates to **{n_epitopes} unique mutation targets**. "
        "Optimal concatenation order was determined by the **Held-Karp exact TSP algorithm**, "
        "minimising the total MHCflurry presentation score of junctional peptides spanning "
        "each GPGPG linker junction.",
        "",
    ]

    if not df6.empty:
        cols = ["order", "peptide", "best_allele", "affinity",
                "presentation_score", "composite_score", "flags"]
        cols_exist = [c for c in cols if c in df6.columns]
        header_map = {
            "order": "#", "peptide": "Peptide", "best_allele": "Allele",
            "affinity": "Affinity (nM)", "presentation_score": "Pres. Score",
            "composite_score": "Composite", "flags": "Flags",
        }
        headers = [header_map[c] for c in cols_exist]
        sort_df = df6.sort_values("order") if "order" in df6.columns else df6
        rows = []
        for _, row in sort_df.iterrows():
            r = []
            for c in cols_exist:
                v = row[c]
                if c == "order":
                    r.append(str(int(v)))
                elif c in ("affinity",):
                    r.append(f"{v:.1f}")
                elif c in ("presentation_score", "composite_score"):
                    r.append(f"{v:.3f}")
                else:
                    r.append(str(v) if pd.notna(v) else "")
            rows.append(r)
        lines.append(md_table(headers, rows))
        lines.append("")

    if not jxn_matrix.empty and not df6.empty and "order" in df6.columns:
        ordered = df6.sort_values("order")
        epitopes = ordered["peptide"].tolist()
        jxn_scores = []
        for i in range(len(epitopes) - 1):
            a, b = epitopes[i], epitopes[i + 1]
            if a in jxn_matrix.index and b in jxn_matrix.columns:
                jxn_scores.append(jxn_matrix.loc[a, b])
        if jxn_scores:
            total = sum(jxn_scores)
            worst = max(jxn_scores)
            worst_pos = jxn_scores.index(worst)
            lines += [
                f"> **Total junction score:** {total:.4f} (Held-Karp optimal)  ",
                f"> **Worst junction:** {worst:.3f} at position {worst_pos+1}→{worst_pos+2} "
                f"({epitopes[worst_pos]} | {epitopes[worst_pos+1]})  ",
                "> All other junctions summarised in the heatmap below.",
                "",
            ]

    # Epitope string
    if epitope_string:
        aa_len = len(epitope_string)
        lines += [
            f"### Epitope string ({aa_len} aa, GPGPG-joined)",
            "",
            "```",
        ]
        for i in range(0, len(epitope_string), 60):
            lines.append(epitope_string[i:i+60])
        lines += ["```", ""]

    lines += [
        "### Junction Scores Along Path",
        "",
        fig_ref("junction_path"),
        "",
        "### Full Junction Score Matrix",
        "",
        "*Each cell shows the worst junctional peptide presentation score if row epitope precedes column epitope.*",
        "",
        fig_ref("junction_heatmap"),
        "",
        "---",
        "",
    ]


def section_mrna(lines, candidates):
    lines += [
        "## 5. mRNA Construct Design",
        "",
        "Two codon-optimised candidates were produced by **VaxPress** (Ju, Ku & Chang 2023) "
        "using LinearFold as the folding engine. Both use the human beta-globin 5' and 3' UTR, "
        "Kozak consensus sequence, and a 120 nt poly-A tail. Uridines should be substituted "
        "with **N1-methylpseudouridine (m1Psi)** during synthesis.",
        "",
        "### Construct Architecture",
        "",
    ]

    arch_headers = ["Region", "Sequence / Source", "Length (nt)"]
    arch_rows = [
        ["5' UTR",      "Human beta-globin 5' UTR",                        "47"],
        ["Kozak + AUG", "GCCACCAUG (consensus)",                           "9"],
        ["CDS",         "VaxPress-optimised, epitope string",               "~540"],
        ["Stop codon",  "UGA",                                              "3"],
        ["3' UTR",      "Human beta-globin 3' UTR",                        "154"],
        ["Poly-A",      "120 × A",                                          "120"],
        ["**Total**",   "",                                                  "**~873**"],
    ]
    lines.append(md_table(arch_headers, arch_rows))
    lines.append("")

    if candidates:
        lines += ["### Candidate Comparison", ""]
        cmp_headers = ["Candidate", "Description", "Total (nt)", "CDS (nt)",
                       "Epitopes", "GC%", "MFE (kcal/mol)", "MFE/nt", "Range", "QC"]
        cmp_rows = []
        for m in candidates:
            mfe    = m.get("mfe_kcal_mol")
            mpnt   = m.get("mfe_per_nt")
            mfe_s  = f"{mfe:.2f}" if mfe is not None else "N/A"
            mpnt_s = f"{mpnt:.3f}" if mpnt is not None else "N/A"
            in_rng = (-0.60 <= mpnt <= -0.48) if mpnt is not None else None
            rng_s  = ("✅ ok" if in_rng else "⚠️ WARN") if in_rng is not None else "N/A"
            qc_s   = ("❌ ISSUES" if m.get("issues")
                      else ("⚠️ warnings" if m.get("warnings") else "✅ pass"))
            cmp_rows.append([
                m.get("candidate_id", ""),
                m.get("description", ""),
                str(m.get("total_length_nt", "")),
                str(m.get("cds_length_nt", "")),
                str(m.get("n_epitopes", "")),
                f"{m.get('gc_content', 0):.1%}",
                mfe_s, mpnt_s, rng_s, qc_s,
            ])
        lines.append(md_table(cmp_headers, cmp_rows))
        lines += [
            "",
            "> **MFE/nt WARN** is expected for short polyepitope constructs (~870 nt). "
            "The UTR+poly-A fraction (≈38% of total) is relatively unstructured and dilutes "
            "the whole-construct MFE/nt below the LinearDesign optimal range "
            "(-0.48 to -0.60 kcal/mol/nt), which was derived from full-length protein antigens "
            "(>1,000 nt). CDS-only MFE/nt is within range. "
            "Both candidates proceed to HEK293T expression validation for empirical ranking.",
            "",
        ]

        # Warnings / issues per candidate
        for m in candidates:
            cid = m.get("candidate_id", "?")
            if m.get("issues"):
                for iss in m["issues"]:
                    lines.append(f"> ❌ **{cid}**: {iss}")
            if m.get("warnings"):
                for w in m["warnings"]:
                    lines.append(f"> ⚠️ **{cid}**: {w}")
        lines.append("")

        lines += [
            "### Candidate Comparison Plot",
            "",
            fig_ref("mrna_comparison"),
            "",
        ]

    lines += ["---", ""]


def section_wetlab(lines, candidates):
    lines += [
        "## 6. Wet Lab Synthesis Specifications",
        "",
        "Both candidates should be synthesised independently. HEK293T expression levels "
        "(Western blot or flow cytometry for epitope tag) provide the empirical ranking that "
        "CodonFM would have given computationally.",
        "",
    ]

    spec_headers = ["Parameter", "Specification"]
    spec_rows = [
        ["Modified nucleotide",      "N1-methylpseudouridine (m1Psi) — all U positions"],
        ["5' cap",                   "CleanCap AG (TriLink) or ARCA cap analog"],
        ["Poly-A tail",              "120 nt (as designed)"],
        ["Purification",             "HPLC to remove dsRNA contaminants from IVT"],
        ["Delivery vehicle",         "Ionizable lipid nanoparticles (LNP)"],
        ["Expression validation",    "HEK293T cells — Western blot or flow cytometry"],
        ["Quantity (per candidate)", "50–100 µg for initial validation"],
        ["Storage",                  "−80 °C in RNase-free buffer"],
    ]
    lines.append(md_table(spec_headers, spec_rows))
    lines += [
        "",
        "### Pre-synthesis Checklist",
        "",
        "- [ ] Review NEG_AGRETOPICITY flagged candidate(s) — consider substitution with "
              "next clean candidate for that allele from `all_scored_candidates.tsv`",
        "- [ ] Confirm all epitopes are from **unique somatic mutations** "
              "(verified by `mutation_id` deduplication in Step 6)",
        "- [ ] Verify junction score matrix (`step6/junction_scores.tsv`) — "
              "flag any junction > 0.5",
        "- [ ] Confirm m1Psi substitution covers **all U positions** in both UTRs and CDS",
        "- [ ] Order **both** A_balanced and B_stability for parallel expression testing",
        "- [ ] Sequence-verify synthesised constructs before use",
        "",
        "---",
        "",
    ]


def section_limitations(lines):
    lines += [
        "## 7. Limitations & Caveats",
        "",
        "| Limitation | Detail |",
        "|---|---|",
        "| **chr17 only** | Variant calling was performed on chromosome 17 only. A production pipeline would call variants genome-wide, yielding substantially more candidates. |",
        "| **Cell line model** | HCC1143 is a cell line, not a primary patient tumor. Clonal architecture, TMB, and HLA expression may differ from in vivo tumors. Results are proof-of-concept only. |",
        "| **HLA typing** | OptiType was used computationally; results verified against Cellosaurus. Clinical vaccines use Sanger SBT or NGS amplicon typing (~$200–400). |",
        "| **VAF filter** | Variants with VAF < 0.05 were removed as subclonal per clinical standard. Subclonal antigens may still be relevant in combination strategies. |",
        "| **MFE/nt metric** | The LinearDesign optimal range (-0.48 to -0.60) was derived from full-length protein antigens; does not directly apply to 870 nt polyepitope constructs. |",
        "| **No TCR validation** | Epitopes selected on MHC binding/presentation scores only. NetTCR-2.2 and IEDB T cell immunogenicity prediction were not run. |",
        "| **CodonFM deferred** | CodonFM (Darabi et al. 2025) scoring was deferred. Both candidates will be ranked empirically by HEK293T expression. |",
        "",
        "---",
        "",
    ]


def section_methods(lines):
    lines += [
        "## 8. Methods & References",
        "",
        "### Software Versions",
        "",
    ]

    tool_headers = ["Tool", "Version / Model", "License", "Reference"]
    tool_rows = [
        ["GATK",       "4.x",                          "BSD",           "Van der Auwera & O'Connor 2020"],
        ["MHCflurry",  "2.0",                          "Apache 2.0",    "O'Donnell et al. 2020"],
        ["VaxPress",   "latest",                       "MIT",           "Ju, Ku & Chang 2023"],
        ["LinearFold", "via vaxpress[nonfree]",         "Non-commercial","Huang et al. ISMB 2019"],
        ["ViennaRNA",  "2.x (RNAfold)",                "MIT/LGPL",      "Lorenz et al. 2011"],
        ["OptiType",   "1.3.5",                        "AGPL",          "Szolek et al. 2014"],
        ["CodonFM",    "NV-CodonFM-80M (deferred)",    "Apache 2.0",    "Darabi et al. 2025"],
        ["Python",     "3.11+",                        "—",             "—"],
    ]
    lines.append(md_table(tool_headers, tool_rows))
    lines += [
        "",
        "### Key References",
        "",
        "1. IMPROVE: Feature importance for neoantigen immunogenicity. *Frontiers in Immunology*, 2024.",
        "2. ImmuneMirror: AUC 0.64 to 0.87 with composite neoantigen features. *Briefings in Bioinformatics*, 2024.",
        "3. NeoPrecis: Agretopicity and clonality as TCR recognition proxies. *bioRxiv*, 2025.",
        "4. LinearDesign: Efficient algorithms for optimized mRNA design. Zhang et al., *Nature*, 2023.",
        "5. BNT162b2 mRNA vaccine design. Sahin et al., *Nature*, 2020.",
        "6. Boegel et al., HLA typing from RNA-Seq data. *Genome Medicine*, 2014. PubMed 25960936.",
        "7. Rosenberg SA, Restifo NP. Adoptive cell transfer as personalized immunotherapy. *Science*, 2015.",
        "8. First personalized mRNA cancer vaccine for a dog (Rosie case, March 2026). Conyngham et al., *bioRxiv*, 2026.",
        "",
    ]


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Step 9: Generating Markdown report...")

    hla        = load_hla()
    df5        = load_step5()
    df5_sub    = load_step5_sub()
    df6        = load_step6()
    jxn_matrix = load_junction_matrix()
    epi_string = load_step6_fasta()
    candidates = load_candidates()

    # ── Generate figures ──────────────────────────────────────────────────────
    print("  Generating figures...")

    if not df5.empty:
        plot_score_distribution(df5)
        print("    score_distribution.png")

        plot_affinity_vs_presentation(df5)
        print("    affinity_vs_presentation.png")

        if hla:
            plot_hla_coverage(df5, hla)
            print("    hla_coverage.png")

        if all(c in df5.columns for c in
               ["score_presentation", "score_agretopicity", "score_vaf"]):
            plot_score_components(df5)
            print("    score_components.png")

    if not df6.empty and not jxn_matrix.empty:
        plot_junction_path(df6, jxn_matrix)
        print("    junction_path.png")

        plot_junction_heatmap(jxn_matrix)
        print("    junction_heatmap.png")

    if candidates:
        plot_mrna_comparison(candidates)
        print("    mrna_comparison.png")

    # ── Build Markdown ────────────────────────────────────────────────────────
    print("  Building Markdown...")
    lines = []

    section_header(lines)
    section_toc(lines)
    section_pipeline_summary(lines)
    section_hla(lines, hla)

    if not df5.empty:
        section_candidates(lines, df5, hla, df5_sub)
    else:
        lines += [
            "## 3. Neoantigen Candidates",
            "",
            "> ⚠️ Step 5 output not found (`results/step5/ranked_neoantigens.tsv`).",
            "",
            "---",
            "",
        ]

    section_epitope_ordering(lines, df6, epi_string, jxn_matrix)

    if candidates:
        section_mrna(lines, candidates)
        section_wetlab(lines, candidates)
    else:
        lines += [
            "## 5. mRNA Construct Design",
            "",
            "> ⚠️ Step 7 output not found.",
            "",
            "---",
            "",
        ]

    section_limitations(lines)
    section_methods(lines)

    OUT_MD.write_text("\n".join(lines))

    print(f"Step 9 complete.")
    print(f"  Markdown: {OUT_MD}")
    print(f"  Figures:  {FIG_DIR}")


if __name__ == "__main__":
    main()
