import os
import re
from setup import agent

BRAIDS_DIR = '/Users/mtwatson/Library/CloudStorage/Box-Box/Projects/AI agent for breeding/breeding data/BRAIDs'
OUTPUTS_DIR = "outputs"

BRAID_EXTENSIONS = {".yaml", ".yml"}

def file_stem_from_path(file_path: str) -> str:
    """Derive a clean file stem from an input filename."""
    stem = os.path.splitext(os.path.basename(file_path))[0]
    # Replace spaces and special characters with underscores
    clean = re.sub(r"[^\w]+", "_", stem).strip("_")
    return clean


def run_on_braid_yaml(braid_program_path: str, output_path: str):
    """Run the agent on a single BRAID YAML file."""
    print(f"\n{'='*60}")
    print(f"BRAID:  {braid_program_path}")
    print(f"Output: {output_path}")
    print(f"{'='*60}")

    task = (
        "Convert the provided BRAID breeding program abstraction YAML into an AlphaSimPy Jupyter notebook (.ipynb).\n\n"
        "The BRAID path is in additional_args as 'braid_program_path'. The notebook output path is in additional_args as 'output_path'.\n\n"
        "Steps:\n"
        "1. Load the BRAID abstraction from braid_program_path.\n"
        "2. Use alphasimpy_tutorials examples as references for notebook structure and format.\n"
        "3. Build a tutorial-style AlphaSimPy notebook with markdown and code cells.\n"
        "4. Save it with save_alphasimpy_notebook(cells=..., output_path=output_path).\n"
        "5. Return the notebook path via final_answer(output_path)."
    )

    try:
        result = agent.run(
            task,
            additional_args={
                "braid_program_path": braid_program_path,
                "output_path": output_path,
            },
        )
        print(f"Done: {result}")
        return result
    except Exception as e:
        print(f"ERROR on {braid_program_path}: {e}")
        return None


def main():
    os.makedirs(OUTPUTS_DIR, exist_ok=True)

    braid_files = sorted(
        f for f in os.listdir(BRAIDS_DIR)
        if os.path.splitext(f)[1].lower() in BRAID_EXTENSIONS
    )

    if not braid_files:
        print(f"No BRAID YAML files found in {BRAIDS_DIR}")
        return

    print(f"Found {len(braid_files)} BRAID file(s) in {BRAIDS_DIR}")

    results = {}
    for filename in braid_files:
        braid_program_path = os.path.join(BRAIDS_DIR, filename)
        file_stem = file_stem_from_path(filename)
        output_path = os.path.join(OUTPUTS_DIR, f"{file_stem}.ipynb")
        results[filename] = run_on_braid_yaml(braid_program_path, output_path)

    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    for filename, result in results.items():
        status = "OK" if result else "FAILED"
        print(f"  [{status}] {filename}")


if __name__ == "__main__":
    main()
