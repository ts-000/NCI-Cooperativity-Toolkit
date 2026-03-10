from Bio.PDB import PDBParser
import itertools
import math
import gzip
import os

def distance(atom1, atom2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(atom1.coord, atom2.coord)))

# Load structure
parser = PDBParser(QUIET=True)
data_dir = os.path.join(os.path.expanduser("~"), "Downloads", "data")
if not os.path.isdir(data_dir):
    raise FileNotFoundError(f"Data directory not found: {data_dir}")

structure = None
for fname in os.listdir(data_dir):
    if fname.endswith(".pdb.gz"):
        with gzip.open(os.path.join(data_dir, fname), "rt") as handle:
            structure = parser.get_structure(os.path.splitext(fname)[0], handle)
        break
if structure is None:
    raise FileNotFoundError(f"No .pdb.gz files found in {data_dir}")

# Aromatic residues to check (π–π interactions)
aromatic_residues = ["PHE", "TYR", "TRP"]

atoms = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_resname() in aromatic_residues:
                for atom in residue:
                    atoms.append((residue.get_resname(), atom))

# Check distances between aromatic atoms
results = []
for (res1, atom1), (res2, atom2) in itertools.combinations(atoms, 2):
    d = distance(atom1, atom2)
    if d < 5.0:  # threshold for π–π stacking
        results.append((res1, atom1.get_name(), res2, atom2.get_name(), round(d, 2)))

# Save results
os.makedirs("results", exist_ok=True)
with open(os.path.join("results", "interactions.csv"), "w") as f:
    f.write("Residue1,Atom1,Residue2,Atom2,Distance\n")
    for r in results:
        f.write(",".join(map(str, r)) + "\n")

print("Found", len(results), "potential π–π interactions.")