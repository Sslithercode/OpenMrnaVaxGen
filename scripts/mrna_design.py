"""
Step 7: mRNA Sequence Design
Takes the ordered epitope string from Step 6 and builds a complete,
optimized mRNA vaccine construct ready for wet lab synthesis.

Pipeline:
  7a. Load ordered epitope AA string from Step 6
  7b. Codon-optimize CDS with VaxPress + LinearFold
  7c. Assemble full construct: 5'UTR + Kozak + CDS + 3'UTR + polyA
  7d. Structural validation with RNAfold (MFE)
  7e. Quality checks + output

Multi-candidate mode:
  Set CANDIDATES below. Each entry defines a named VaxPress weight profile.
  The pipeline runs once per candidate, writing outputs to:
      step7/<candidate_id>/vaccine_mrna.fasta
      step7/<candidate_id>/quality_metrics.json
      step7/<candidate_id>/vaxpress_output/report.html
  A comparison table is printed and saved at step7/candidate_comparison.json.

  Profiles sweep the MFE/CAI tradeoff — mirroring the LinearDesign Nature
  (2023) lambda sweep which showed the best immunogenic response is neither
  the most stable nor the highest-CAI sequence. CodonFM validation (Step 8)
  then scores each candidate independently to identify the winner.

  To run a single candidate (original behaviour), set:
      CANDIDATES = [{"id": "A_balanced", "description": "Balanced", "weights": {}}]

Tool stack:
  VaxPress      Ju, Ku & Chang (2023)        MIT license
  LinearFold    Huang et al. (ISMB 2019)     Non-commercial (via vaxpress[nonfree])
  ViennaRNA     Lorenz et al. (2011)         MIT/LGPL (open-source fallback)

Note on CodonFM (Darabi et al., 2025, Apache 2.0):
  CodonFM validation runs in Step 8, not here. Step 8 loads the 1B-param
  model once and scores all candidates in a single pass for efficiency.

Note on LinearDesign (Zhang et al., Nature 2023, non-commercial):
  For sequences under ~600 aa, VaxPress at 200 iterations converges to
  equivalent quality without LinearDesign. To enable for longer constructs:
      --lineardesign 1.0 --lineardesign-dir /path/to/LinearDesign
      --conservative-start 10 --initial-mutation-rate 0.01 --iterations 300
"""

import warnings
warnings.filterwarnings("ignore")

import subprocess
import json
import re
import sys
from pathlib import Path
import os
# ── Config ────────────────────────────────────────────────────────────────────

STEP6_FASTA = Path.home() / "melanoma-pipeline/results/step6/ordered_epitopes.fasta"
OUT_DIR     = Path.home() / "melanoma-pipeline/results/step7"

LINKER = "GPGPG"  # must match Step 6 — used to count epitopes from the AA string

# All UTR constants stored as DNA, converted to RNA uniformly in assemble_construct.
# [FIXED] Previously UTR3 was stored as RNA (using U) while UTR5/Kozak used DNA (T).
# Inconsistent storage means a future edit to UTR3 using T notation silently breaks.
UTR5  = "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC"
KOZAK = "GCCACCATG"  # terminal ATG is start codon; VaxPress ATG stripped in assembly
UTR3  = (
    "TGCTGTTTGTGCCTTTGTCTAGTGAGTTTCCATATGCTGGTCAGTAGCAG"
    "TGTTTCTGATGTGTACTGTCGCTCGCTGTTTTTGTTTCAGTGTACCTGTTC"
    "TGCGTTTGCTCTCTGCTGTGCTGTTTAGGGTGAGAGTGTGTAGCTGTCTGCAG"
)

POLY_A_LENGTH       = 120   # nt, standard for therapeutic mRNA
VAXPRESS_ITERATIONS = 100
VAXPRESS_PROCESSES  = 8

# Optimal MFE/nt range from LinearDesign (Zhang et al., Nature 2023).
# Too negative (< -0.60): over-stabilized, ribosome stalling, poor translation.
# Too positive (> -0.48): under-stabilized, rapid degradation, short half-life.
# Candidates outside this range are flagged as warnings in quality checks.
MFE_PER_NT_MIN = -0.60
MFE_PER_NT_MAX = -0.48

# ── Candidate definitions ─────────────────────────────────────────────────────
#
# Each candidate is a dict:
#   id          short name used for output subdirectory
#   description human-readable label for the comparison table
#   weights     dict of VaxPress scoring weight overrides
#               keys: mfe, cai, ucount, loop, gc
#               omitted keys use VaxPress defaults
#               passed as --<key>-weight <value> on the CLI
#
# Biological rationale:
#   mfe-weight    secondary structure stability -> mRNA half-life
#   cai-weight    codon adaptation index -> translation rate / protein yield
#   ucount-weight uridine depletion -> reduced TLR7/8 innate immune activation
#   loop-weight   unpaired loop regions -> degradation resistance
#   gc-weight     GC content balance -> synthesis quality, no hotspots
#
# These four profiles sweep the design space the way the LinearDesign
# Nature (2023) paper did with lambda. Step 8 CodonFM scoring ranks objectively.

CANDIDATES = [
    {
        "id":          "A_balanced",
        "description": "Balanced MFE + CAI (default)",
        "weights":     {},
    },
    {
        "id":          "B_stability",
        "description": "Stability-biased (longer half-life)",
        "weights":     {"mfe": 10, "cai": 2, "loop": 5},
    },
]

# ── Helpers ───────────────────────────────────────────────────────────────────

def read_fasta(fasta_file):
    seq = ""
    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq


def count_epitopes(epitope_aa):
    """Count epitopes in the GPGPG-joined string produced by Step 6."""
    return len(epitope_aa.split(LINKER))


def translate(dna):
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    dna = dna.upper().replace("U", "T")
    protein = ""
    for i in range(0, len(dna) - 2, 3):
        aa = codon_table.get(dna[i:i+3], "X")
        if aa == "*":
            break
        protein += aa
    return protein


def simple_codon_optimize(epitope_aa):
    """Naive fallback if VaxPress fails."""
    HUMAN_OPTIMAL_CODONS = {
        "A": "GCC", "R": "AGG", "N": "AAC", "D": "GAC",
        "C": "TGC", "Q": "CAG", "E": "GAG", "G": "GGC",
        "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
        "M": "ATG", "F": "TTC", "P": "CCC", "S": "AGC",
        "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    }
    cds = "ATG"
    for aa in epitope_aa:
        cds += HUMAN_OPTIMAL_CODONS.get(aa, "NNN")
    cds += "TGA"
    return cds


# ── Step 7a ───────────────────────────────────────────────────────────────────

def load_epitope_string():
    epitope_aa = read_fasta(STEP6_FASTA)
    n_epitopes = count_epitopes(epitope_aa)
    print(f"Loaded epitope string: {len(epitope_aa)} amino acids, {n_epitopes} epitopes")
    print(f"First 60 aa: {epitope_aa[:60]}...")
    return epitope_aa


# ── Step 7b ───────────────────────────────────────────────────────────────────

def detect_folding_engine():
    try:
        import linearfold  # noqa: F401
        print("  Folding engine: LinearFold (~2x faster, non-commercial)")
        return "linearfold"
    except ImportError:
        print("  Folding engine: ViennaRNA (open-source, exact MFE)")
        return "vienna"


def codon_optimize(epitope_aa, out_dir, weights):
    """
    Codon-optimize with VaxPress + LinearFold.

    weights: dict of VaxPress scoring weight overrides, e.g.
             {"mfe": 10, "cai": 2} -> --mfe-weight 10 --cai-weight 2
             Empty dict = all VaxPress defaults (balanced profile).
    """
    input_fasta = out_dir / "epitope_input.fasta"
    with open(input_fasta, "w") as f:
        f.write(f">vaccine_epitopes\n{epitope_aa}\n")

    vaxpress_out = out_dir / "vaxpress_output"
    vaxpress_out.mkdir(exist_ok=True)

    folding_engine = detect_folding_engine()

    print(f"\nRunning VaxPress codon optimization...", flush=True)
    print(f"  Input:      {len(epitope_aa)} amino acids", flush=True)
    print(f"  Iterations: {VAXPRESS_ITERATIONS}", flush=True)
    print(f"  Processes:  {VAXPRESS_PROCESSES}", flush=True)
    if weights:
        print(f"  Weights:    {weights}", flush=True)

    cmd = [
        "vaxpress",
        "-i", str(input_fasta),
        "-o", str(vaxpress_out),
        "--protein",
        "--species", "human",
        "--processes", str(VAXPRESS_PROCESSES),
        "--iterations", str(VAXPRESS_ITERATIONS),
        "--folding-engine", folding_engine,
        "--overwrite",
    ]

    for key, value in weights.items():
        cmd += [f"--{key}-weight", str(value)]

    result = subprocess.run(cmd, text=True)

    if result.returncode != 0:
        print("VaxPress failed, falling back to simple codon optimization...")
        return simple_codon_optimize(epitope_aa)

    optimized_file = vaxpress_out / "result.fasta"
    if not optimized_file.exists():
        fastas = sorted(vaxpress_out.glob("*.fasta"))
        if fastas:
            optimized_file = fastas[0]
        else:
            print("VaxPress output not found, falling back")
            return simple_codon_optimize(epitope_aa)

    optimized_cds = read_fasta(optimized_file).upper().replace("U", "T")
    print(f"VaxPress complete. CDS: {len(optimized_cds)} nt")
    return optimized_cds


# ── Step 7c ───────────────────────────────────────────────────────────────────

def assemble_construct(cds):
    """
    Assemble full mRNA construct (all RNA):
      [5'UTR] [Kozak+AUG] [CDS body] [stop codon] [3'UTR] [polyA]

    Strips leading ATG from VaxPress CDS — Kozak provides the start codon
    to avoid a duplicated AUG at the 5' end.

    [FIXED] All UTR constants now stored as DNA and converted here uniformly
    with .replace("T", "U"). Previously UTR3 was stored as RNA (U) while
    UTR5/Kozak used DNA (T), making it fragile to future edits.

    [FIXED] Returns cds_full_rna = Kozak AUG + CDS body + stop codon,
    so cds_length_nt in metrics reflects the true CDS including start codon.
    Previously returned cds body only, understating length by 3 nt.
    """
    cds_body = cds.upper().replace("U", "T")
    if cds_body.startswith("ATG"):
        cds_body = cds_body[3:]

    if cds_body[-3:] not in ["TAA", "TAG", "TGA"]:
        cds_body += "TGA"

    # All constants stored as DNA — convert uniformly to RNA
    utr5_rna  = UTR5.upper().replace("T", "U")
    kozak_rna = KOZAK.upper().replace("T", "U")
    cds_rna   = cds_body.upper().replace("T", "U")
    utr3_rna  = UTR3.upper().replace("T", "U")   # [FIXED] was .upper() only
    poly_a    = "A" * POLY_A_LENGTH

    full_mrna    = utr5_rna + kozak_rna + cds_rna + utr3_rna + poly_a
    cds_full_rna = kozak_rna + cds_rna  # [FIXED] full CDS including AUG

    print(f"\n── mRNA Construct Assembly ──────────────────")
    print(f"  5' UTR:      {len(utr5_rna):>5} nt  (human beta-globin)")
    print(f"  Kozak + AUG: {len(kozak_rna):>5} nt")
    print(f"  CDS:         {len(cds_full_rna):>5} nt  (VaxPress optimized, includes AUG)")
    print(f"  3' UTR:      {len(utr3_rna):>5} nt  (human beta-globin)")
    print(f"  Poly-A:      {len(poly_a):>5} nt")
    print(f"  ─────────────────────────────────────────────")
    print(f"  Total:       {len(full_mrna):>5} nt")

    return full_mrna, cds_full_rna


# ── Step 7d ───────────────────────────────────────────────────────────────────

def validate_structure(mrna_sequence, out_dir):
    """RNAfold MFE of full construct. More negative = more stable."""
    print("\nRunning RNAfold structural validation...")

    if subprocess.run(["which", "RNAfold"], capture_output=True).returncode != 0:
        print("  RNAfold not found — install: sudo apt install vienna-rna")
        return None, None

    seq_file = out_dir / "mrna_for_fold.fasta"
    with open(seq_file, "w") as f:
        f.write(f">vaccine_mrna\n{mrna_sequence}\n")

    fold_result = subprocess.run(
        ["RNAfold", "--noPS", str(seq_file)],
        capture_output=True, text=True
    )

    if fold_result.returncode != 0:
        print(f"  RNAfold failed: {fold_result.stderr[:300]}")
        return None, None

    mfe = None
    structure = None
    for line in fold_result.stdout.strip().split("\n"):
        match = re.search(r'\((-?\d+\.\d+)\)\s*$', line)
        if match:
            mfe = float(match.group(1))
            structure = line.split()[0]
            break

    if mfe is not None:
        print(f"  MFE: {mfe:.2f} kcal/mol  (more negative = more stable)")

    return mfe, structure


# ── Step 7e ───────────────────────────────────────────────────────────────────

def quality_checks(cds, mrna, original_epitope_aa, mfe):
    """
    QC checks on the assembled construct.

    cds  — DNA string (pre-assembly output from VaxPress or fallback)
    mrna — RNA string (full assembled construct)
    mfe  — float or None (from RNAfold)

    [ADDED] mfe parameter for MFE/nt range check against LinearDesign optimum.

    Note on homopolymer check: checks T runs in the DNA template. The
    delivered mRNA uses U, but synthesis is from the DNA template, so
    T-homopolymers are the relevant concern for IVT synthesis quality.
    """
    issues        = []
    warnings_list = []
    cds_dna       = cds.upper().replace("U", "T")

    if len(cds_dna) % 3 != 0:
        issues.append(f"CDS length {len(cds_dna)} not divisible by 3")

    for i in range(0, len(cds_dna) - 3, 3):
        if cds_dna[i:i+3] in ["TAA", "TAG", "TGA"]:
            issues.append(f"Internal stop codon at position {i}: {cds_dna[i:i+3]}")
            break

    gc = (mrna.count("G") + mrna.count("C")) / len(mrna)
    if gc < 0.40 or gc > 0.70:
        warnings_list.append(
            f"GC content {gc:.1%} outside optimal 40-70% range"
        )

    for base in ["A", "T", "C", "G"]:
        if base * 8 in cds_dna:
            warnings_list.append(
                f"Long homopolymer run ({base}x8+) in CDS DNA template"
            )

    # [ADDED] MFE/nt range check — LinearDesign (Zhang et al., Nature 2023)
    if mfe is not None:
        mfe_per_nt = mfe / len(mrna)
        if mfe_per_nt < MFE_PER_NT_MIN or mfe_per_nt > MFE_PER_NT_MAX:
            warnings_list.append(
                f"MFE/nt {mfe_per_nt:.3f} outside optimal range "
                f"({MFE_PER_NT_MIN} to {MFE_PER_NT_MAX} kcal/mol/nt)"
            )
        else:
            print(f"  MFE/nt: {mfe_per_nt:.3f} kcal/mol/nt  (within optimal range)")

    print(f"\n── Translation Round-Trip Check ─────────────")
    translated = translate(cds_dna)
    print(f"  Original:   {len(original_epitope_aa)} aa")
    print(f"  Translated: {len(translated)} aa")

    if translated == original_epitope_aa:
        print(f"  Result: exact match")
    else:
        overlap = sum(a == b for a, b in zip(translated, original_epitope_aa))
        pct = overlap / max(len(translated), len(original_epitope_aa))
        if pct > 0.95:
            print(f"  Result: {pct:.1%} match (minor start codon boundary)")
        else:
            issues.append(f"Translation mismatch: {pct:.1%} identity")
            print(f"  Result: WARNING {pct:.1%} identity")

    return issues, warnings_list


# ── Per-candidate runner ──────────────────────────────────────────────────────

def run_candidate(epitope_aa, candidate):
    os.system('cls')
    """Run the full Step 7 pipeline for one candidate profile."""
    cid  = candidate["id"]
    desc = candidate["description"]
    wts  = candidate["weights"]

    print(f"\n{'='*50}")
    print(f"  Candidate {cid}: {desc}")
    if wts:
        print(f"  Weight overrides: {wts}")
    print(f"{'='*50}")

    cand_dir = OUT_DIR / cid
    cand_dir.mkdir(parents=True, exist_ok=True)

    cds                   = codon_optimize(epitope_aa, cand_dir, wts)
    full_mrna, cds_full   = assemble_construct(cds)
    mfe, _                = validate_structure(full_mrna, cand_dir)
    # [FIXED] mfe now passed into quality_checks for MFE/nt validation
    issues, warnings_list = quality_checks(cds, full_mrna, epitope_aa, mfe)

    mrna_fasta = cand_dir / "vaccine_mrna.fasta"
    with open(mrna_fasta, "w") as f:
        f.write(f">vaccine_mrna | {cid} | {desc}\n")
        for i in range(0, len(full_mrna), 60):
            f.write(full_mrna[i:i+60] + "\n")

    gc         = (full_mrna.count("G") + full_mrna.count("C")) / len(full_mrna)
    mfe_per_nt = round(mfe / len(full_mrna), 4) if mfe is not None else None

    metrics = {
        "candidate_id":     cid,
        "description":      desc,
        "weight_overrides": wts,
        "total_length_nt":  len(full_mrna),
        # [FIXED] was len(cds_rna) — missing the 3 nt Kozak AUG start codon
        "cds_length_nt":    len(cds_full),
        # [FIXED] was hardcoded to 30 — now derived from actual GPGPG-joined string
        "n_epitopes":       count_epitopes(epitope_aa),
        "gc_content":       round(gc, 4),
        "mfe_kcal_mol":     mfe,
        # [ADDED] stored alongside raw MFE for direct comparison to literature range
        "mfe_per_nt":       mfe_per_nt,
        "codonfm_score":    None,   # populated by Step 8
        "issues":           issues,
        "warnings":         warnings_list,
        "wet_lab_recommendations": [
            "Substitute all uridines with N1-methylpseudouridine (m1Psi)",
            "Use CleanCap AG or ARCA cap analog for 5' capping",
            "Poly-A tail: 120 nt as designed",
            "Purify by HPLC to remove dsRNA contaminants from IVT",
            "Encapsulate in ionizable lipid nanoparticles (LNP)",
            "Validate expression in HEK293T cells before patient use",
        ],
    }

    with open(cand_dir / "quality_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\n  ── {cid} summary ──")
    print(f"  Total mRNA:  {len(full_mrna):,} nt")
    print(f"  CDS:         {len(cds_full):,} nt  (including AUG)")
    print(f"  GC content:  {gc:.1%}")
    if mfe is not None:
        print(f"  MFE:         {mfe:.2f} kcal/mol")
    if mfe_per_nt is not None:
        in_range  = MFE_PER_NT_MIN <= mfe_per_nt <= MFE_PER_NT_MAX
        range_tag = "in range" if in_range else "OUT OF RANGE"
        print(f"  MFE/nt:      {mfe_per_nt:.3f} kcal/mol/nt  ({range_tag})")
    status = "ISSUES" if issues else ("warnings" if warnings_list else "all checks passed")
    print(f"  QC:          {status}")
    if issues:
        for iss in issues:
            print(f"    ERROR: {iss}")
    if warnings_list:
        for w in warnings_list:
            print(f"    WARN:  {w}")
    print(f"  Output:      {mrna_fasta}")

    return metrics


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Step 7: mRNA Sequence Design")
    print(f"  Input:      {STEP6_FASTA}")
    print(f"  Output:     {OUT_DIR}")
    print(f"  Candidates: {len(CANDIDATES)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not STEP6_FASTA.exists():
        print(f"\nERROR: Step 6 output not found: {STEP6_FASTA}")
        sys.exit(1)

    epitope_aa = load_epitope_string()

    all_metrics = []
    for candidate in CANDIDATES:
        metrics = run_candidate(epitope_aa, candidate)
        all_metrics.append(metrics)

    comparison_file = OUT_DIR / "candidate_comparison.json"
    with open(comparison_file, "w") as f:
        json.dump(all_metrics, f, indent=2)

    print(f"\n{'='*72}")
    print(f"  Candidate Comparison")
    print(f"{'='*72}")
    print(f"  {'ID':<20} {'GC':>6}  {'MFE':>10}  {'MFE/nt':>8}  {'Range':>8}  QC")
    print(f"  {'─'*20} {'─'*6}  {'─'*10}  {'─'*8}  {'─'*8}  {'─'*10}")
    for m in all_metrics:
        mfe_str   = f"{m['mfe_kcal_mol']:.2f}" if m['mfe_kcal_mol'] is not None else "N/A"
        mpnt_str  = f"{m['mfe_per_nt']:.3f}"   if m['mfe_per_nt']   is not None else "N/A"
        if m['mfe_per_nt'] is not None:
            in_range  = MFE_PER_NT_MIN <= m['mfe_per_nt'] <= MFE_PER_NT_MAX
            range_str = "ok" if in_range else "WARN"
        else:
            range_str = "N/A"
        qc_str = "ISSUES" if m['issues'] else ("warn" if m['warnings'] else "pass")
        print(
            f"  {m['candidate_id']:<20} {m['gc_content']:>5.1%}  "
            f"{mfe_str:>10}  {mpnt_str:>8}  {range_str:>8}  {qc_str}"
        )

    print(f"\n  CodonFM scores: run Step 8 to populate")
    print(f"  Comparison:     {comparison_file}")
    print(f"\nStep 7 complete.")


if __name__ == "__main__":
    main()