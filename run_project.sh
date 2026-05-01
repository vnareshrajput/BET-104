#!/bin/bash
set -uo pipefail

if [ -z "${1:-}" ]; then
    echo "Usage: ./run_project.sh /path/to/pdb_folder [cores]"
    exit 1
fi

PDB_DIR=$(realpath "$1")
CORES=${2:-$(sysctl -n hw.ncpu)}

mkdir -p extracted stride_files results/angles results/plots

snakemake \
    --snakefile snakemake_workflow.smk \
    --cores "$CORES" \
    --rerun-incomplete \
    --keep-going \
    --config pdb_dir="$PDB_DIR"

python scripts/make_plot.py results/angles results/plots/final_plot.png
echo "Done."
