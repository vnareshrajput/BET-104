#!/usr/bin/env python3
"""
plot_angles.py
Usage: python3 plot_angles.py <angle_dir> <output_path>

Dependencies (install via pip):
    pip install numpy matplotlib scipy
"""

import sys
import os
import glob
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend — safe on macOS without a display
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import gaussian_kde


SIZE_ORDER = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]

COLORS = {
    "Tiny":         "#F5F0E8",
    "Small":        "#E8C97A",
    "Intermediate": "#D4845A",
    "Large":        "#B84040",
    "Bulky":        "#7A0C0C",
}


def load_angles(angle_dir):
    data = defaultdict(list)
    valid_count = 0

    # Use Path.glob for reliable cross-platform behaviour on macOS
    for fpath in Path(angle_dir).glob("*.txt"):
        if fpath.stat().st_size == 0:
            continue

        file_had_data = False

        with fpath.open() as f:
            next(f, None)  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                try:
                    angle = float(parts[2])
                    size = parts[4]
                    if size not in SIZE_ORDER:
                        continue
                    data[size].append(angle)
                    file_had_data = True
                except Exception:
                    continue

        if file_had_data:
            valid_count += 1

    return data, valid_count


def plot(data, valid_count, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    total = sum(len(v) for v in data.values())

    x = np.linspace(-180, 180, 1000)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Background styling
    fig.patch.set_facecolor("#8C8C8C")
    ax.set_facecolor("#8C8C8C")

    # vertical grid lines
    for xv in np.arange(-150, 181, 25):
        ax.axvline(xv, color="#6688CC", linestyle="dotted", linewidth=0.7, alpha=0.8)

    # KDE curves
    for size in SIZE_ORDER:
        angles = np.array(data[size])
        if len(angles) < 10:
            continue

        kde = gaussian_kde(angles, bw_method=0.15)
        ax.plot(x, kde(x), color=COLORS[size], linewidth=2.5, label=size)

    # axes styling
    ax.set_xlim(-180, 180)
    ax.set_xlabel("Angle between adjacent C-α -> Centroid vectors [°]", fontsize=13)
    ax.set_ylabel("Norm. Freq. [A.U.]", fontsize=13)

    ax.set_title(
        f"Tripeptide (XRX) in Helix (n = {total})",
        fontsize=14,
        pad=10
    )

    ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.tick_params(axis="both", labelsize=11)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    ax.legend(
        frameon=True,
        fontsize=12,
        facecolor="white",
        edgecolor="white",
        loc="upper left"
    )

    plt.tight_layout()
    plt.savefig(str(output_path), dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()

    print(f"Saved: {output_path}")


def main():
    angle_dir = sys.argv[1]
    output_path = sys.argv[2]

    data, valid_count = load_angles(angle_dir)

    total = sum(len(v) for v in data.values())

    if total == 0:
        print("No data found.")
        sys.exit(1)

    for size in SIZE_ORDER:
        print(f"{size}: {len(data[size])}")

    plot(data, valid_count, output_path)


if __name__ == "__main__":
    main()