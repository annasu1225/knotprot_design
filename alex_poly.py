'''
Example usage:
returns a dict of knot probabilities if Closure.TWO_POINTS is used
for deterministic restults, use Closure.MASS_CENTER
'''

from topoly import alexander, Closure
import argparse as ap

def alex_poly(pdb, tries, max_cross):
    return alexander(pdb, Closure.TWO_POINTS, tries=tries, max_cross=max_cross)

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("pdb_file", help="Path to the PDB file")
    parser.add_argument("--tries", type=int, default=100, help="Number of closure attempts")
    parser.add_argument("--max_cross", type=int, default=25, help="Maximum crossings to consider")
    args = parser.parse_args()
    knot_prob = alex_poly(args.pdb_file, tries=args.tries, max_cross=args.max_cross)
    print(knot_prob)

# Example usage:
# python alex_poly.py path/to/pdbfile.pdb --tries 200 --max_cross 10