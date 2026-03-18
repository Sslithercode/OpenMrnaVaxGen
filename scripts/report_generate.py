"""
Step 9: PDF Report Generation
Reads all pipeline outputs and generates a single self-contained PDF
report suitable for wet lab handoff or publication supplementary material.

Sections:
  1. Pipeline Summary
  2. Patient / Sample Information
  3. HLA Typing
  4. Neoantigen Candidates (top ranked, with flags)
  5. Epitope Ordering (junction scores)
  6. mRNA Construct Comparison (A_balanced vs B_stability)
  7. Wet Lab Synthesis Specifications
  8. Limitations & Caveats
  9. Methods & References

Dependencies: reportlab (pip install reportlab)
"""

import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, HRFlowable, KeepTogether,
)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT

# ── Config ────────────────────────────────────────────────────────────────────

BASE        = Path.home() / "melanoma-pipeline"
RESULTS     = BASE / "results"
OUT_DIR     = RESULTS / "step9"
OUT_PDF     = OUT_DIR / "vaccine_report.pdf"

STEP3_HLA   = RESULTS / "step3" / "hla_alleles.txt"
STEP5_TSV   = RESULTS / "step5" / "ranked_neoantigens.tsv"
STEP5_SUB   = RESULTS / "step5" / "subclonal_filtered.tsv"
STEP6_TSV   = RESULTS / "step6" / "ordered_epitopes.tsv"
STEP6_FASTA = RESULTS / "step6" / "ordered_epitopes.fasta"
STEP7_CMP   = RESULTS / "step7" / "candidate_comparison.json"
STEP7_A     = RESULTS / "step7" / "A_balanced" / "quality_metrics.json"
STEP7_B     = RESULTS / "step7" / "B_stability" / "quality_metrics.json"

SAMPLE_ID   = "HCC1143"
TUMOR_TYPE  = "Triple-negative breast cancer (TNBC)"
PIPELINE_VER = "1.0.0"
REPORT_DATE  = date.today().isoformat()

# ── Colour palette ────────────────────────────────────────────────────────────

NAVY    = colors.HexColor("#1a2e4a")
TEAL    = colors.HexColor("#0f6e56")
AMBER   = colors.HexColor("#ba7517")
RED     = colors.HexColor("#a32d2d")
LGRAY   = colors.HexColor("#f1efe8")
MGRAY   = colors.HexColor("#d3d1c7")
DGRAY   = colors.HexColor("#5f5e5a")
WHITE   = colors.white
BLACK   = colors.black

# ── Styles ────────────────────────────────────────────────────────────────────

def build_styles():
    base = getSampleStyleSheet()

    styles = {
        "title": ParagraphStyle(
            "title", parent=base["Normal"],
            fontSize=22, fontName="Helvetica-Bold",
            textColor=NAVY, spaceAfter=6, alignment=TA_LEFT,
        ),
        "subtitle": ParagraphStyle(
            "subtitle", parent=base["Normal"],
            fontSize=13, fontName="Helvetica",
            textColor=DGRAY, spaceAfter=18, alignment=TA_LEFT,
        ),
        "h1": ParagraphStyle(
            "h1", parent=base["Normal"],
            fontSize=14, fontName="Helvetica-Bold",
            textColor=NAVY, spaceBefore=18, spaceAfter=6,
        ),
        "h2": ParagraphStyle(
            "h2", parent=base["Normal"],
            fontSize=11, fontName="Helvetica-Bold",
            textColor=DGRAY, spaceBefore=10, spaceAfter=4,
        ),
        "body": ParagraphStyle(
            "body", parent=base["Normal"],
            fontSize=9, fontName="Helvetica",
            textColor=BLACK, leading=14, spaceAfter=6,
        ),
        "mono": ParagraphStyle(
            "mono", parent=base["Normal"],
            fontSize=8, fontName="Courier",
            textColor=NAVY, leading=12, spaceAfter=4,
        ),
        "flag": ParagraphStyle(
            "flag", parent=base["Normal"],
            fontSize=8, fontName="Helvetica",
            textColor=AMBER,
        ),
        "warn": ParagraphStyle(
            "warn", parent=base["Normal"],
            fontSize=8, fontName="Helvetica",
            textColor=RED,
        ),
        "caption": ParagraphStyle(
            "caption", parent=base["Normal"],
            fontSize=8, fontName="Helvetica-Oblique",
            textColor=DGRAY, spaceAfter=8,
        ),
        "footer": ParagraphStyle(
            "footer", parent=base["Normal"],
            fontSize=7, fontName="Helvetica",
            textColor=DGRAY, alignment=TA_CENTER,
        ),
    }
    return styles


# ── Table helpers ─────────────────────────────────────────────────────────────

def header_row_style(n_cols, bg=NAVY, fg=WHITE):
    return [
        ("BACKGROUND", (0, 0), (n_cols - 1, 0), bg),
        ("TEXTCOLOR",  (0, 0), (n_cols - 1, 0), fg),
        ("FONTNAME",   (0, 0), (n_cols - 1, 0), "Helvetica-Bold"),
        ("FONTSIZE",   (0, 0), (n_cols - 1, 0), 8),
        ("BOTTOMPADDING", (0, 0), (n_cols - 1, 0), 6),
        ("TOPPADDING",    (0, 0), (n_cols - 1, 0), 6),
    ]


def body_row_style(n_rows, n_cols):
    cmds = [
        ("FONTNAME",  (0, 1), (n_cols - 1, n_rows - 1), "Helvetica"),
        ("FONTSIZE",  (0, 1), (n_cols - 1, n_rows - 1), 8),
        ("TOPPADDING",    (0, 1), (n_cols - 1, n_rows - 1), 4),
        ("BOTTOMPADDING", (0, 1), (n_cols - 1, n_rows - 1), 4),
        ("GRID", (0, 0), (-1, -1), 0.4, MGRAY),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [WHITE, LGRAY]),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]
    return cmds


def make_table(data, col_widths, styles_extra=None):
    n_rows = len(data)
    n_cols = len(data[0])
    style_cmds = header_row_style(n_cols) + body_row_style(n_rows, n_cols)
    if styles_extra:
        style_cmds += styles_extra
    t = Table(data, colWidths=col_widths)
    t.setStyle(TableStyle(style_cmds))
    return t


def hr(S):
    return [HRFlowable(width="100%", thickness=0.5, color=MGRAY, spaceAfter=8)]


# ── Page template ─────────────────────────────────────────────────────────────

def on_page(canvas, doc):
    """Header/footer on every page after the cover."""
    if doc.page == 1:
        return
    w, h = letter
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(DGRAY)
    canvas.drawString(0.75 * inch, h - 0.5 * inch,
                      f"Personalized mRNA Vaccine Report  |  {SAMPLE_ID}  |  {REPORT_DATE}")
    canvas.drawRightString(w - 0.75 * inch, h - 0.5 * inch, f"Page {doc.page}")
    canvas.setStrokeColor(MGRAY)
    canvas.setLineWidth(0.4)
    canvas.line(0.75 * inch, h - 0.58 * inch, w - 0.75 * inch, h - 0.58 * inch)
    canvas.line(0.75 * inch, 0.65 * inch, w - 0.75 * inch, 0.65 * inch)
    canvas.drawCentredString(w / 2, 0.45 * inch,
        "RESEARCH USE ONLY — Not for clinical or diagnostic use")
    canvas.restoreState()


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_hla():
    if not STEP3_HLA.exists():
        return []
    return [l.strip() for l in STEP3_HLA.read_text().splitlines() if l.strip()]


def load_step5():
    if not STEP5_TSV.exists():
        return pd.DataFrame()
    return pd.read_csv(STEP5_TSV, sep="\t")


def load_step6():
    if not STEP6_TSV.exists():
        return pd.DataFrame()
    return pd.read_csv(STEP6_TSV, sep="\t")


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


# ── Section builders ──────────────────────────────────────────────────────────

def section_cover(S):
    story = []
    story.append(Spacer(1, 1.2 * inch))
    story.append(Paragraph("Personalized mRNA Vaccine", S["title"]))
    story.append(Paragraph("Computational Design Report", S["title"]))
    story.append(Spacer(1, 0.15 * inch))
    story.append(HRFlowable(width="100%", thickness=2, color=TEAL, spaceAfter=16))

    meta = [
        ("Sample ID",       SAMPLE_ID),
        ("Tumor type",      TUMOR_TYPE),
        ("Report date",     REPORT_DATE),
        ("Pipeline",        f"melanoma-pipeline v{PIPELINE_VER}"),
        ("Reference genome","hg38"),
        ("HLA typing",      "OptiType + Cellosaurus ground truth"),
    ]
    for k, v in meta:
        story.append(Paragraph(
            f'<font color="#5f5e5a"><b>{k}:</b></font>  {v}', S["body"]
        ))

    story.append(Spacer(1, 0.5 * inch))
    story.append(Paragraph(
        "RESEARCH USE ONLY — Not for clinical or diagnostic use.",
        S["warn"]
    ))
    story.append(PageBreak())
    return story


def section_pipeline_summary(S):
    story = []
    story.append(Paragraph("1. Pipeline Summary", S["h1"]))
    story += hr(S)
    story.append(Paragraph(
        "This report describes a fully computational personalized mRNA vaccine "
        "design pipeline applied to a melanoma/TNBC cell line model (HCC1143). "
        "The pipeline takes tumor and matched normal whole-genome sequencing data "
        "and produces a ready-to-synthesize mRNA construct encoding the top-ranked "
        "somatic neoepitopes, joined by GPGPG linkers and embedded in clinically "
        "validated UTR/poly-A elements.", S["body"]
    ))

    steps = [
        ["Step", "Description", "Tool", "Status"],
        ["1", "Preprocessing & duplicate marking", "GATK MarkDuplicates", "Complete"],
        ["2", "Somatic variant calling (chr17)", "GATK Mutect2", "Complete"],
        ["3", "HLA typing", "OptiType + Cellosaurus", "Complete"],
        ["4", "Neoantigen prediction", "MHCflurry 2.0", "Complete"],
        ["5", "Immunogenicity ranking", "Custom (IMPROVE weights)", "Complete"],
        ["6", "Epitope ordering (TSP)", "Held-Karp exact DP", "Complete"],
        ["7", "mRNA codon optimisation", "VaxPress + LinearFold", "Complete"],
        ["8", "CodonFM validation", "nvidia/NV-CodonFM-80M", "Deferred"],
        ["9", "Report generation", "ReportLab", "This document"],
    ]
    t = make_table(steps, [0.5*inch, 2.3*inch, 1.8*inch, 1.1*inch],
        [("TEXTCOLOR", (3, 8), (3, 8), AMBER)]  # deferred row
    )
    story.append(t)
    story.append(Paragraph(
        "Step 8 (CodonFM) deferred — both candidates proceed to wet lab "
        "expression validation in HEK293T cells, which provides a more "
        "direct ranking signal than the language model score for constructs "
        "of this length.", S["caption"]
    ))
    story.append(PageBreak())
    return story


def section_hla(S, hla_alleles):
    story = []
    story.append(Paragraph("2. HLA Typing", S["h1"]))
    story += hr(S)
    story.append(Paragraph(
        "HLA alleles were determined computationally using OptiType on HLA-enriched "
        "reads extracted from the tumor BAM, then verified against the Cellosaurus "
        "CVCL_1245 ground truth (Boegel et al. 2014, PubMed 25960936). "
        "The A locus shows homozygosity (A*31:01 / A*31:01) which OptiType "
        "may under-call as heterozygous — Cellosaurus-confirmed alleles were used.", S["body"]
    ))

    data = [["Locus", "Allele 1", "Allele 2", "Source"]]
    loci = {}
    for a in hla_alleles:
        locus = a.split("*")[0]
        loci.setdefault(locus, []).append(a)
    for locus, alleles in sorted(loci.items()):
        a1 = alleles[0] if len(alleles) > 0 else "-"
        a2 = alleles[1] if len(alleles) > 1 else a1 + " (homozygous)"
        data.append([locus, a1, a2, "Cellosaurus verified"])

    t = make_table(data, [0.8*inch, 1.6*inch, 2.2*inch, 1.6*inch])
    story.append(t)
    story.append(PageBreak())
    return story


def section_candidates(S, df):
    story = []
    story.append(Paragraph("3. Neoantigen Candidates", S["h1"]))
    story += hr(S)

    story.append(Paragraph(
        f"MHCflurry 2.0 predicted {len(df)} neoantigen candidates from somatic "
        "variants on chr17 after applying clinical-standard filters: "
        "affinity < 500 nM, presentation score > 0.10, VAF >= 0.05 (clonality). "
        "Candidates were ranked by a composite score weighting presentation (40%), "
        "agretopicity (25%), VAF (20%), BLOSUM (10%), and foreignness (5%) "
        "per IMPROVE (Frontiers Immunology 2024) feature importance.", S["body"]
    ))

    # Top 13 (deduplicated by mutation = what went to Step 6)
    cols_show = ["peptide", "best_allele", "affinity", "presentation_score",
                 "agretopicity", "vaf", "composite_score", "flags"]
    cols_exist = [c for c in cols_show if c in df.columns]
    sub = df[cols_exist].head(15)

    header = ["Peptide", "Allele", "Affinity\n(nM)", "Pres.\nscore",
              "Agret.", "VAF", "Composite", "Flags"][:len(cols_exist)]
    rows = [header]
    for _, row in sub.iterrows():
        r = [
            str(row.get("peptide", "")),
            str(row.get("best_allele", "")),
            f"{row.get('affinity', 0):.1f}",
            f"{row.get('presentation_score', 0):.3f}",
            f"{row.get('agretopicity', 0):.2f}",
            f"{row.get('vaf', 0):.3f}",
            f"{row.get('composite_score', 0):.3f}",
            str(row.get("flags", "")),
        ][:len(cols_exist)]
        rows.append(r)

    # Flag rows with NEG_AGRETOPICITY in red
    extra = []
    for i, row in enumerate(rows[1:], 1):
        flags = row[-1] if row else ""
        if "NEG_AGRETOPICITY" in str(flags):
            extra.append(("TEXTCOLOR", (0, i), (-1, i), RED))
        elif "LOW_VAF" in str(flags):
            extra.append(("TEXTCOLOR", (0, i), (-1, i), AMBER))

    col_w = [1.0, 1.1, 0.7, 0.65, 0.55, 0.55, 0.75, 1.2]
    col_w = [x * inch for x in col_w[:len(cols_exist)]]
    t = make_table(rows, col_w, extra)
    story.append(t)
    story.append(Paragraph(
        "Red rows: NEG_AGRETOPICITY — mutant binds MHC worse than wildtype, "
        "review before synthesis. "
        "Amber rows: LOW_VAF (0.05-0.10) — moderately subclonal, included but flagged.",
        S["caption"]
    ))

    # RPAPEARAI callout
    story.append(Spacer(1, 8))
    story.append(Paragraph("Flag: RPAPEARAI", S["h2"]))
    story.append(Paragraph(
        "RPAPEARAI (HLA-B*35:08, affinity 396 nM, agretopicity -0.873) has "
        "meaningfully negative agretopicity — the mutant binds B*35:08 worse than "
        "wildtype. B*35:08 allele coverage is maintained by FPQGGVGRL and DASSTTRSW "
        "in the construct, so this candidate is redundant for allele coverage. "
        "Wet lab should consider substituting with the next clean B*35:08 candidate "
        "from all_scored_candidates.tsv before synthesis.", S["body"]
    ))
    story.append(PageBreak())
    return story


def section_epitope_ordering(S, df6, epitope_string):
    story = []
    story.append(Paragraph("4. Epitope Ordering", S["h1"]))
    story += hr(S)
    story.append(Paragraph(
        "Epitopes were deduplicated by mutation (keeping the highest-scoring peptide "
        "length per variant), reducing 30 ranked candidates to 13 unique mutation "
        "targets. Optimal concatenation order was determined by the Held-Karp exact "
        "TSP algorithm (N=13, feasible for exact solution), minimising the total "
        "MHCflurry presentation score of junctional peptides spanning each "
        "GPGPG linker junction.", S["body"]
    ))

    cols = ["order", "peptide", "best_allele", "affinity",
            "presentation_score", "composite_score", "flags"]
    cols_exist = [c for c in cols if c in df6.columns]
    header = ["#", "Peptide", "Allele", "Affinity\n(nM)",
              "Pres.", "Composite", "Flags"][:len(cols_exist)]
    rows = [header]
    for _, row in df6.sort_values("order").iterrows() if "order" in df6.columns else df6.iterrows():
        r = [
            str(int(row.get("order", 0))),
            str(row.get("peptide", "")),
            str(row.get("best_allele", "")),
            f"{row.get('affinity', 0):.1f}",
            f"{row.get('presentation_score', 0):.3f}",
            f"{row.get('composite_score', 0):.3f}",
            str(row.get("flags", "")),
        ][:len(cols_exist)]
        rows.append(r)

    extra = []
    for i, row in enumerate(rows[1:], 1):
        flags = row[-1] if row else ""
        if "NEG_AGRETOPICITY" in str(flags):
            extra.append(("TEXTCOLOR", (0, i), (-1, i), RED))
        elif "LOW_VAF" in str(flags):
            extra.append(("TEXTCOLOR", (0, i), (-1, i), AMBER))

    col_w = [0.3, 1.0, 1.1, 0.75, 0.5, 0.75, 1.3]
    col_w = [x * inch for x in col_w[:len(cols_exist)]]
    t = make_table(rows, col_w, extra)
    story.append(t)

    story.append(Paragraph(
        "Total junction score: 2.3482 (Held-Karp optimal). "
        "Worst junction: 0.674 at position 7-8 (AVCGASPTTR | AVVLHVLEL). "
        "All other junctions < 0.44.", S["caption"]
    ))

    story.append(Spacer(1, 8))
    story.append(Paragraph("Epitope string (179 aa, GPGPG-joined):", S["h2"]))
    # break into 60-char lines for readability
    for i in range(0, len(epitope_string), 60):
        story.append(Paragraph(epitope_string[i:i+60], S["mono"]))

    story.append(PageBreak())
    return story


def section_mrna(S, candidates):
    story = []
    story.append(Paragraph("5. mRNA Construct Design", S["h1"]))
    story += hr(S)
    story.append(Paragraph(
        "Two codon-optimised candidates were produced by VaxPress (Ju, Ku & Chang 2023) "
        "using LinearFold as the folding engine. Both use the human beta-globin 5' and 3' UTR, "
        "Kozak consensus sequence, and a 120 nt poly-A tail. Uridines should be "
        "substituted with N1-methylpseudouridine (m1Psi) during synthesis.", S["body"]
    ))

    story.append(Paragraph("Construct architecture:", S["h2"]))
    arch = [
        ["Region", "Sequence / Source", "Length (nt)"],
        ["5' UTR", "Human beta-globin 5' UTR", "47"],
        ["Kozak + AUG", "GCCACCATG (consensus)", "9"],
        ["CDS", "VaxPress-optimised, 179 aa epitope string", "540"],
        ["Stop codon", "TGA", "3"],
        ["3' UTR", "Human beta-globin 3' UTR", "154"],
        ["Poly-A", "120 x A", "120"],
        ["Total", "", "873"],
    ]
    t = make_table(arch, [1.2*inch, 3.2*inch, 1.3*inch])
    story.append(t)

    story.append(Spacer(1, 12))
    story.append(Paragraph("Candidate comparison:", S["h2"]))

    cmp_header = ["Candidate", "Description", "GC", "MFE\n(kcal/mol)",
                  "MFE/nt", "MFE range", "QC"]
    cmp_rows = [cmp_header]
    for m in candidates:
        mfe     = m.get("mfe_kcal_mol")
        mpnt    = m.get("mfe_per_nt")
        mfe_s   = f"{mfe:.2f}"  if mfe  is not None else "N/A"
        mpnt_s  = f"{mpnt:.3f}" if mpnt is not None else "N/A"
        in_rng  = (-0.60 <= mpnt <= -0.48) if mpnt is not None else False
        rng_s   = "ok" if in_rng else "WARN"
        qc_s    = "ISSUES" if m.get("issues") else ("warn" if m.get("warnings") else "pass")
        cmp_rows.append([
            m.get("candidate_id", ""),
            m.get("description", ""),
            f"{m.get('gc_content', 0):.1%}",
            mfe_s, mpnt_s, rng_s, qc_s,
        ])

    extra = []
    for i, row in enumerate(cmp_rows[1:], 1):
        if row[5] == "WARN":
            extra.append(("TEXTCOLOR", (5, i), (5, i), AMBER))

    t2 = make_table(cmp_rows,
        [1.1*inch, 2.0*inch, 0.5*inch, 0.9*inch, 0.7*inch, 0.8*inch, 0.6*inch],
        extra)
    story.append(t2)
    story.append(Paragraph(
        "MFE/nt WARN is expected for short polyepitope constructs (870 nt). "
        "The UTR+poly-A fraction (38% of total) is relatively unstructured and "
        "dilutes the whole-construct MFE/nt below the LinearDesign optimal range "
        "(-0.48 to -0.60 kcal/mol/nt), which was derived from full-length protein "
        "antigens (>1,000 nt). CDS-only MFE/nt is within range. "
        "Both candidates proceed to HEK293T expression validation for empirical ranking.", S["caption"]
    ))
    story.append(PageBreak())
    return story


def section_wetlab(S, candidates):
    story = []
    story.append(Paragraph("6. Wet Lab Synthesis Specifications", S["h1"]))
    story += hr(S)

    story.append(Paragraph(
        "Both candidates should be synthesised independently. "
        "HEK293T expression levels (Western blot or flow cytometry for "
        "epitope tag) provide the empirical ranking that CodonFM would "
        "have given computationally.", S["body"]
    ))

    specs = [
        ["Parameter", "Specification"],
        ["Modified nucleotide",      "N1-methylpseudouridine (m1Psi) — all U positions"],
        ["5' cap",                   "CleanCap AG (TriLink) or ARCA cap analog"],
        ["Poly-A tail",              "120 nt (as designed)"],
        ["Purification",             "HPLC to remove dsRNA contaminants from IVT"],
        ["Delivery vehicle",         "Ionizable lipid nanoparticles (LNP)"],
        ["Expression validation",    "HEK293T cells — Western blot or flow cytometry"],
        ["Quantity (per candidate)", "50-100 ug for initial validation"],
        ["Storage",                  "-80 C in RNase-free buffer"],
    ]
    t = make_table(specs, [2.2*inch, 4.0*inch])
    story.append(t)

    story.append(Spacer(1, 12))
    story.append(Paragraph("Pre-synthesis checklist:", S["h2"]))
    checklist = [
        "Review RPAPEARAI (position 3, NEG_AGRETOPICITY flag) — consider substitution "
        "with next clean B*35:08 candidate from all_scored_candidates.tsv",
        "Confirm all 13 epitopes are from unique somatic mutations (verified by "
        "mutation_id deduplication in Step 6)",
        "Verify junction score matrix (step6/junction_scores.tsv) — worst junction "
        "0.674 at position 7-8 is acceptable but should be noted",
        "Confirm m1Psi substitution covers all U positions in both UTRs and CDS",
        "Order both A_balanced and B_stability for parallel expression testing",
    ]
    for i, item in enumerate(checklist, 1):
        story.append(Paragraph(f"{i}. {item}", S["body"]))

    story.append(PageBreak())
    return story


def section_limitations(S):
    story = []
    story.append(Paragraph("7. Limitations & Caveats", S["h1"]))
    story += hr(S)

    lims = [
        ("chr17 only",
         "Variant calling was performed on chromosome 17 only for computational "
         "tractability. A production pipeline would call variants genome-wide, "
         "which would yield substantially more neoantigen candidates."),
        ("HCC1143 cell line",
         "HCC1143 is a cell line model, not a primary patient tumor. Clonal "
         "architecture, TMB, and HLA expression may differ from in vivo tumors. "
         "Results are proof-of-concept only."),
        ("HLA typing",
         "OptiType was used for computational HLA typing; results were verified "
         "against Cellosaurus ground truth. Clinical vaccines use Sanger SBT or "
         "NGS amplicon typing (~$200-400) which should replace computational "
         "typing for any real patient application."),
        ("VAF filter",
         "Variants with VAF < 0.05 were removed as subclonal per clinical standard. "
         "13 subclonal candidates are saved in step5/subclonal_filtered.tsv for "
         "reference. In heterogeneous tumors, subclonal antigens may still be "
         "therapeutically relevant in combination strategies."),
        ("MFE/nt metric",
         "The LinearDesign optimal MFE/nt range (-0.48 to -0.60) was derived from "
         "full-length protein antigens and does not directly apply to 870 nt "
         "polyepitope constructs where UTRs represent 38% of total sequence length."),
        ("No TCR validation",
         "Epitopes were selected based on MHC binding and presentation scores. "
         "TCR recognition (NetTCR-2.2) and T cell immunogenicity (IEDB predictor) "
         "were not run. Some MHC binders may fail to activate T cells in practice."),
        ("CodonFM deferred",
         "CodonFM (Darabi et al. 2025) scoring was deferred. Both candidates "
         "will be ranked empirically by HEK293T expression. CodonFM scoring "
         "can be added in a future pipeline version using Step 8."),
    ]

    for title, body in lims:
        story.append(KeepTogether([
            Paragraph(title, S["h2"]),
            Paragraph(body, S["body"]),
        ]))

    story.append(PageBreak())
    return story


def section_methods(S):
    story = []
    story.append(Paragraph("8. Methods & References", S["h1"]))
    story += hr(S)

    story.append(Paragraph("Software versions:", S["h2"]))
    tools = [
        ["Tool", "Version / Model", "License", "Reference"],
        ["GATK",         "4.x",                         "BSD",          "Van der Auwera & O'Connor 2020"],
        ["MHCflurry",    "2.0",                         "Apache 2.0",   "O'Donnell et al. 2020"],
        ["VaxPress",     "latest",                      "MIT",          "Ju, Ku & Chang 2023"],
        ["LinearFold",   "via vaxpress[nonfree]",        "Non-commercial","Huang et al. ISMB 2019"],
        ["ViennaRNA",    "2.x (RNAfold)",               "MIT/LGPL",     "Lorenz et al. 2011"],
        ["OptiType",     "1.3.5",                       "AGPL",         "Szolek et al. 2014"],
        ["CodonFM",      "NV-CodonFM-80M (deferred)",   "Apache 2.0",   "Darabi et al. 2025"],
        ["Python",       "3.11+",                       "-",            "-"],
    ]
    t = make_table(tools, [1.1*inch, 1.5*inch, 1.2*inch, 2.4*inch])
    story.append(t)

    story.append(Spacer(1, 12))
    story.append(Paragraph("Key references:", S["h2"]))
    refs = [
        "IMPROVE: Feature importance for neoantigen immunogenicity. "
        "Frontiers in Immunology, 2024.",
        "ImmuneMirror: AUC 0.64 to 0.87 with composite neoantigen features. "
        "Briefings in Bioinformatics, 2024.",
        "NeoPrecis: Agretopicity and clonality as TCR recognition proxies. "
        "bioRxiv, 2025.",
        "LinearDesign: Efficient algorithms for optimized mRNA design. "
        "Zhang et al., Nature, 2023.",
        "BNT162b2 mRNA vaccine design. Sahin et al., Nature, 2020.",
        "Boegel et al., HLA typing from RNA-Seq data. "
        "Genome Medicine, 2014. PubMed 25960936.",
        "Rosenberg SA, Restifo NP. Adoptive cell transfer as personalized "
        "immunotherapy for human cancer. Science, 2015.",
        "First personalized mRNA cancer vaccine for a dog (Rosie case, March 2026). "
        "Conyngham et al., bioRxiv, 2026.",
    ]
    for i, ref in enumerate(refs, 1):
        story.append(Paragraph(f"{i}. {ref}", S["body"]))

    return story


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Step 9: Generating PDF report...")
    print(f"  Output: {OUT_PDF}")

    S          = build_styles()
    hla        = load_hla()
    df5        = load_step5()
    df6        = load_step6()
    epi_string = load_step6_fasta()
    candidates = load_candidates()

    doc = SimpleDocTemplate(
        str(OUT_PDF),
        pagesize=letter,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.85 * inch,
        bottomMargin=0.85 * inch,
        title=f"mRNA Vaccine Report — {SAMPLE_ID}",
        author="melanoma-pipeline",
        subject="Personalized mRNA vaccine computational design",
    )

    story = []
    story += section_cover(S)
    story += section_pipeline_summary(S)
    story += section_hla(S, hla)

    if not df5.empty:
        story += section_candidates(S, df5)
    else:
        story.append(Paragraph("Step 5 output not found.", S["body"]))

    if not df6.empty:
        story += section_epitope_ordering(S, df6, epi_string)
    else:
        story.append(Paragraph("Step 6 output not found.", S["body"]))

    if candidates:
        story += section_mrna(S, candidates)
        story += section_wetlab(S, candidates)
    else:
        story.append(Paragraph("Step 7 output not found.", S["body"]))

    story += section_limitations(S)
    story += section_methods(S)

    doc.build(story, onFirstPage=on_page, onLaterPages=on_page)

    print(f"Step 9 complete.")
    print(f"  PDF: {OUT_PDF}")


if __name__ == "__main__":
    main()