import random
import argparse as ap
from topoly import alexander, Closure
import os
from datetime import datetime
from pathlib import Path

from esm.sdk.api import ESMProtein
from esm.sdk.experimental import GuidedDecodingScoringFunction
from esm.utils.structure.protein_chain import ProteinChain
from esm.sdk.experimental import (
    ConstraintType,
    ESM3GuidedDecodingWithConstraints,
    GenerationConstraint
)


# --- Access Model ---
# Option 1: Use forge API to access the model
# model = client(
#     model="esm3-medium-2024-08", url="https://forge.evolutionaryscale.ai", token='6Zk4FIijMlhj5iMQqWUU2q'
# )

# Option 2: use local installation of esm
from esm.models.esm3 import ESM3
model = ESM3.from_pretrained().to("cuda")

# --- Define Scoring Functions ---
class PTMScoringFunction(GuidedDecodingScoringFunction):
    def __call__(self, protein: ESMProtein) -> float:
        assert protein.ptm is not None, "Protein must have pTM scores to be scored"
        return float(protein.ptm)
    
class KnotScoringFunction(GuidedDecodingScoringFunction):
    def __init__(self, target_knot, out_dir):
             super().__init__()
             self.target_knot = target_knot
             # output directory for PDBs
             self.out_dir = out_dir

    def __call__(self, protein: ESMProtein) -> float:
        score = self.knot_probability(protein)
        return score

    def knot_probability(self, protein: ESMProtein) -> float:
        print("current working directory:", os.getcwd())
        pdb_f = f"{self.target_knot}.pdb"
        pdb_path = os.path.join(self.out_dir, pdb_f)
        _ = protein.to_pdb(pdb_path)
        try:
            result = alexander(pdb_path, Closure.TWO_POINTS, tries=200, max_cross=5)
            p_knot = result.get(self.target_knot, 0.0)
            print(f"******Protein seq for knot {self.target_knot} with probability {p_knot}: {protein.sequence}")
            return p_knot
        except Exception as e:
            print(f"Alexander polynomial calculation failed: {e}")
            return 0.0

def get_masked_sequence(wildtype, masking_percentage):
    wildtype_pdb = wildtype[0]
    wildtype_chain = wildtype[1]
    wildtype = ESMProtein.from_protein_chain(
            ProteinChain.from_rcsb(wildtype_pdb, chain_id=wildtype_chain)
        )
    maskable_indices = list(range(len(wildtype.sequence)))
    num_to_mask = int(len(maskable_indices) * masking_percentage)
    indices_to_mask = random.sample(maskable_indices, num_to_mask)
    refinement_template_list = list(wildtype.sequence)

    for i in indices_to_mask: 
        refinement_template_list[i] = '_'
    masked_seq = "".join(refinement_template_list)
    print(f"Masked sequence: {masked_seq}")
    return masked_seq

def run_guided_generate(out_dir, seq_len, knot_type, wildtype, num_decoding_steps, num_samples_per_step, masking_percentage):

    if wildtype is not None:
        # Start from a wildtype sequence
        masked_seq = get_masked_sequence(wildtype, masking_percentage)
        starting_protein = ESMProtein(sequence=masked_seq)
    else:
        # Start from a fully masked protein
        starting_protein = ESMProtein(sequence="_" * seq_len)

    # Constrain generation to have pTM > 0.75
    ptm_constraint = GenerationConstraint(
        scoring_function=PTMScoringFunction(),
        constraint_type=ConstraintType.GREATER_EQUAL,
        value=0.75,
    )

    scoring_fn = KnotScoringFunction(knot_type, out_dir)
    knot_guided_decoding = ESM3GuidedDecodingWithConstraints(
        client=model, 
        scoring_function=scoring_fn,
        constraints=[ptm_constraint],
        damping=1.0,  # Damping factor for the MMDM algorithm
        learning_rate=10.0,  # Learning rate for the MMDM algorithm
        )

    # Call guided_generate
    generated_protein = knot_guided_decoding.guided_generate(
        protein=starting_protein,
        num_decoding_steps=num_decoding_steps,
        num_samples_per_step=num_samples_per_step
    )

    # Save the trajectory visualization of the guided decoding
    png_path = os.path.join(out_dir, "trajectory.png")
    import matplotlib.pyplot as plt
    knot_guided_decoding.visualize_latest_trajectory()
    fig = plt.gcf()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return generated_protein

def main():
    parser = ap.ArgumentParser()
    parser.add_argument("--out_dir", type=str, help="Output directory to save generated PDBs")
    parser.add_argument("--seq_length", default=256, type=int, help="Sequence Length for the generated protein")
    parser.add_argument("--knot_type", default="0_1", type=str, help="The knot type of the generated protein")
    parser.add_argument("--wildtype", nargs=2, metavar=("PDB_ID", "CHAIN_ID"), default=None, help="Provide PDB ID and chain (e.g. --wildtype 1uak A) to initiate the guided generation")
    parser.add_argument("--num_decoding_steps", default=64, type=int, help="Number of steps in the guided generation")
    parser.add_argument("--num_samples_per_step", default=10, type=int, help="Number of sample proteins generated in each step")
    parser.add_argument("--masking_percentage", default=0.4, type=float, help="Percentage of the sequence to mask")

    args = parser.parse_args()

    # Build a default out_dir if none provided: runs/<knot_type>_YYYYmmdd_HHMMSS
    if not args.out_dir:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_dir = Path("runs") / f"{args.knot_type}_{ts}"
        out_dir = str(default_dir)
    else:
        out_dir = args.out_dir

    # Ensure output directory exists
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out_dir}")

    # Call guided generation
    run_guided_generate(
        out_dir,
        args.seq_length,
        args.knot_type,
        args.wildtype,
        args.num_decoding_steps,
        args.num_samples_per_step,
        args.masking_percentage
    )

if __name__=='__main__':
    main()

# Example usage:
# python guided_generation.py --knot_type 3_1 --wildtype 1cmx A 
# python guided_generation.py --seq_length 300 --knot_type 3_1