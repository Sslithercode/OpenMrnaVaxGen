"""
Step 4: Neoantigen Prediction — orchestrator
Runs Step 4a (peptide extraction) then Step 4b (MHC binding) as two
separate Python processes to avoid the htslib/TensorFlow double-free crash.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from paths import step_dir

SCRIPTS_DIR = Path(__file__).parent


def main():
    # Step 4a: varcode/pyensembl — no TF in this process
    result = subprocess.run([sys.executable, str(SCRIPTS_DIR / "neoantigen_peptides.py")])
    if result.returncode != 0:
        print("Step 4a (peptide extraction) failed.")
        sys.exit(1)

    # Step 4b: MHCflurry — no htslib in this process
    result = subprocess.run([sys.executable, str(SCRIPTS_DIR / "neoantigen_binding.py")])
    if result.returncode != 0:
        print("Step 4b (binding prediction) failed.")
        sys.exit(1)

    print("\nStep 4 complete.")


if __name__ == "__main__":
    main()
