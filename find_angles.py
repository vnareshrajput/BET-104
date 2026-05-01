#!/usr/bin/env python3
"""
compute_angles_from_tsv.py
Usage: python3 compute_angles_from_tsv.py <tsv_file> <pdb_file> <output_file>

Dependencies (install via pip):
    pip install biopython numpy
"""

import sys
import os
import warnings
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

BACKBONE = {"N", "CA", "C", "O", "OXT"}

SIZE_CATEGORY = {
    "GLY": "Tiny",
    "ALA": "Tiny",
    "SER": "Small",
    "CYS": "Small",
    "THR": "Small",
    "VAL": "Small",
    "LEU": "Intermediate",
    "ILE": "Intermediate",
    "ASN": "Intermediate",
    "ASP": "Intermediate",
    "PRO": "Intermediate",
    "GLN": "Large",
    "GLU": "Large",
    "MET": "Large",
    "HIS": "Large",
    "LYS": "Large",
    "PHE": "Bulky",
    "TYR": "Bulky",
    "TRP": "Bulky",
    "ARG": "Bulky",
}


def side_chain_centroid(residue):
    coords = [
        a.get_vector().get_array()
        for a in residue.get_atoms()
        if a.get_name() not in BACKBONE and not a.get_name().startswith("H")
    ]
    if not coords:
        return residue["CA"].get_vector().get_array()
    return np.mean(coords, axis=0)


def signed_angle(v1, v2, axis):
    axis = np.array(axis, dtype=float)
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)

    if np.linalg.norm(axis) < 1e-8:
        return None

    axis /= np.linalg.norm(axis)
    v1 = v1 - np.dot(v1, axis) * axis   # project out axis component
    v2 = v2 - np.dot(v2, axis) * axis

    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-8 or n2 < 1e-8:
        return None

    v1 /= n1
    v2 /= n2

    return np.degrees(np.arctan2(np.dot(np.cross(v1, v2), axis), np.dot(v1, v2)))


def build_residue_index(structure):
    index = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                key = (chain.id, residue.get_id()[1])
                index[key] = residue
        break
    return index


def main():
    tsv_file = Path(sys.argv[1])
    pdb_file = Path(sys.argv[2])
    output_file = Path(sys.argv[3])

    if tsv_file.stat().st_size == 0:
        output_file.open("w").close()
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", str(pdb_file))
    residue_index = build_residue_index(structure)

    results = []

    with tsv_file.open() as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

            chain = parts[0]
            prev_resnum = int(parts[1])
            prev_resname = parts[2]
            center_resnum = int(parts[4])
            next_resnum = int(parts[7])

            prev_res = residue_index.get((chain, prev_resnum))
            center_res = residue_index.get((chain, center_resnum))
            next_res = residue_index.get((chain, next_resnum))

            if not all([prev_res, center_res, next_res]):
                continue

            try:
                ca_prev = prev_res["CA"].get_vector().get_array()
                ca_center = center_res["CA"].get_vector().get_array()
                centroid_prev = side_chain_centroid(prev_res)
                centroid_center = side_chain_centroid(center_res)
            except KeyError:
                continue

            v1 = ca_prev - centroid_prev
            v2 = ca_center - centroid_center
            axis = ca_center - ca_prev

            angle = signed_angle(v1, v2, axis)
            if angle is None:
                continue

            size = SIZE_CATEGORY.get(prev_resname, "Unknown")
            results.append(f"{chain}\t{center_resnum}\t{angle:.4f}\t{prev_resname}\t{size}")

    with output_file.open("w") as out:
        out.write("chain\tcenter_resnum\tangle\tprev_resname\tsize_category\n")
        for row in results:
            out.write(row + "\n")


if __name__ == "__main__":
    main()