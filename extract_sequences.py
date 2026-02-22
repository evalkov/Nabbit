#!/usr/bin/env python3
"""Extract chain sequences from all PDB files in the same directory. Double-click to run."""

import os, glob

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}

def extract_chains(pdb_path):
    chains = {}
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip()
            if resname not in AA3TO1:
                continue
            chain = line[21].strip() or "_"
            resseq = line[22:27].strip()
            key = (chain, resseq)
            if key in seen:
                continue
            seen.add(key)
            chains.setdefault(chain, []).append(AA3TO1[resname])
    return {c: "".join(residues) for c, residues in chains.items()}

script_dir = os.path.dirname(os.path.abspath(__file__))
pdbs = sorted(glob.glob(os.path.join(script_dir, "*.pdb")))

if not pdbs:
    print("No .pdb files found in", script_dir)
    input("Press Enter to exit...")
    raise SystemExit

out_path = os.path.join(script_dir, "pdb_sequences.tsv")
n = 0
with open(out_path, "w") as out:
    out.write("file\tchain\tlength\tsequence\n")
    for pdb in pdbs:
        fname = os.path.basename(pdb)
        for chain, seq in sorted(extract_chains(pdb).items()):
            out.write(f"{fname}\t{chain}\t{len(seq)}\t{seq}\n")
            n += 1

print(f"Wrote {n} chains from {len(pdbs)} PDB files -> {out_path}")
input("Press Enter to exit...")
