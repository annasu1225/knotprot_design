## KnotProt Design Toolkit

This workspace combines guided sequence generation with knot-aware scoring to design proteins that target specific mathematical knot types. The workflow hinges on two scripts:

- `guided_generation.py` — interfaces with ESM3 to propose sequences conditioned on knot type and pTM constraints.
- `alex_poly.py` — wraps TopoLy's Alexander polynomial routine to score how likely a generated backbone adopts a desired knot.

### Environment Setup

Create and activate the conda environment defined in `esm3_topoly.yml`:

```bash
conda env create -f esm3_topoly.yml
conda activate esm3_topoly
```

The environment ships all required dependencies, including ESM3 and TopoLy.

### Generate Knot-Conditioned Proteins

`guided_generation.py` runs ESM3's guided decoding loop with optional wild-type starting points. Common flags:

- `--knot_type` (e.g., `3_1`, `4_1`, `0_1`)
- `--seq_length` for fully masked starts
- `--wildtype <PDB_ID> <CHAIN_ID>` to mask a known structure
- `--masking_percentage`, `--num_decoding_steps`, `--num_samples_per_step`

Example:

```bash
python guided_generation.py \
	--knot_type 3_1 \
	--seq_length 280 \
	--num_decoding_steps 64 \
	--num_samples_per_step 10
```

Outputs (PDBs, trajectory PNGs, logs) appear under `runs/<knot_type>_<timestamp>/`.

### Score Knot Probabilities

To evaluate any backbone structure (generated or external), call `alex_poly.py` with a PDB path and optional sampling parameters:

```bash
python alex_poly.py generated_examples/3_1_all_atom.pdb --tries 200 --max_cross 20
```

The script prints a dictionary mapping knot identifiers (e.g., `3_1`) to estimated probabilities derived from repeated closure attempts.

### Typical Workflow

1. Generate candidate sequences/backbones with `guided_generation.py` for the desired knot type.
2. (Optional) Run ProteinMPNN inverse folding on promising backbones.
3. Score the resulting PDBs via `alex_poly.py` to confirm topological fidelity.
4. Visualize the PDBs in KnotPlot and run dynamical relaxation. See my other repo [Molecular_Topology](https://github.com/annasu1225/Molecular_Topology) for instructions and codes for this. 

### Acknowledgement
I appreciate the helpful discussions with Prof. Jeffrey Brock, Prof. Mark Gerstein, Prof. Smita Krishnaswamy, Neil Voss, Santanu Antu, and Hiren Madnu.


