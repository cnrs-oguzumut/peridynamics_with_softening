#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def conf_index(path: Path) -> int:
    stem = path.stem
    return int(stem.replace("conf", ""))


def load_conf(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 3:
        raise ValueError(f"{path} must contain position, displacement, strain columns")
    if data.shape[1] == 3:
        damage = np.zeros((data.shape[0], 1))
        data = np.hstack((data, damage))
    return data[:, :4]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot strain profiles from result/conf*.dat files."
    )
    parser.add_argument(
        "--result-dir",
        default="result",
        type=Path,
        help="Directory containing conf*.dat files.",
    )
    parser.add_argument(
        "--steps",
        nargs="*",
        type=int,
        help="Specific configuration numbers to plot, e.g. --steps 1 20 80.",
    )
    parser.add_argument(
        "--displacement",
        action="store_true",
        help="Also plot displacement profiles.",
    )
    parser.add_argument(
        "--damage",
        action="store_true",
        help="Also plot damage profiles.",
    )
    parser.add_argument(
        "--save",
        type=Path,
        help="Save figure to this file instead of only showing it.",
    )
    args = parser.parse_args()

    files = sorted(args.result_dir.glob("conf*.dat"), key=conf_index)
    if args.steps:
        wanted = set(args.steps)
        files = [path for path in files if conf_index(path) in wanted]

    if not files:
        raise SystemExit(f"No conf*.dat files found in {args.result_dir}")

    ncols = 1 + int(args.displacement) + int(args.damage)
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 4), squeeze=False)
    strain_ax = axes[0, 0]
    next_ax = 1
    disp_ax = axes[0, next_ax] if args.displacement else None
    if args.displacement:
        next_ax += 1
    damage_ax = axes[0, next_ax] if args.damage else None

    for path in files:
        data = load_conf(path)
        x = data[:, 0]
        u = data[:, 1]
        strain = data[:, 2]
        damage = data[:, 3]
        label = f"conf {conf_index(path)}"

        strain_ax.plot(x, strain, label=label)
        if disp_ax is not None:
            disp_ax.plot(x, u, label=label)
        if damage_ax is not None:
            damage_ax.plot(x, damage, label=label)

    strain_ax.set_xlabel("Position")
    strain_ax.set_ylabel("Strain")
    strain_ax.set_title("Strain Field")
    strain_ax.grid(True)
    strain_ax.legend()

    if disp_ax is not None:
        disp_ax.set_xlabel("Position")
        disp_ax.set_ylabel("Displacement")
        disp_ax.set_title("Displacement Field")
        disp_ax.grid(True)
        disp_ax.legend()

    if damage_ax is not None:
        damage_ax.set_xlabel("Position")
        damage_ax.set_ylabel("Damage")
        damage_ax.set_title("Damage Field")
        damage_ax.set_ylim(-0.05, 1.05)
        damage_ax.grid(True)
        damage_ax.legend()

    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=200)
    plt.show()


if __name__ == "__main__":
    main()
