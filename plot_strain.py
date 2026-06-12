#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from matplotlib.animation import FFMpegWriter
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


def make_axes(ncols: int) -> tuple[plt.Figure, np.ndarray]:
    return plt.subplots(1, ncols, figsize=(6 * ncols, 4), squeeze=False)


def selected_axes(
    axes: np.ndarray,
    show_displacement: bool,
    show_damage: bool,
    show_strain_zoom: bool,
):
    strain_ax = axes[0, 0]
    next_ax = 1
    zoom_ax = axes[0, next_ax] if show_strain_zoom else None
    if show_strain_zoom:
        next_ax += 1
    disp_ax = axes[0, next_ax] if show_displacement else None
    if show_displacement:
        next_ax += 1
    damage_ax = axes[0, next_ax] if show_damage else None
    return strain_ax, zoom_ax, disp_ax, damage_ax


def style_axes(strain_ax, zoom_ax, disp_ax, damage_ax) -> None:
    strain_ax.set_xlabel("Position")
    strain_ax.set_ylabel("Strain")
    strain_ax.set_title("Strain Field")
    strain_ax.grid(True)

    if zoom_ax is not None:
        zoom_ax.set_xlabel("Position")
        zoom_ax.set_ylabel("Strain")
        zoom_ax.set_title("Strain Field Zoom")
        zoom_ax.grid(True)

    if disp_ax is not None:
        disp_ax.set_xlabel("Position")
        disp_ax.set_ylabel("Displacement")
        disp_ax.set_title("Displacement Field")
        disp_ax.grid(True)

    if damage_ax is not None:
        damage_ax.set_xlabel("Position")
        damage_ax.set_ylabel("Damage")
        damage_ax.set_title("Damage Field")
        damage_ax.set_ylim(-0.05, 1.05)
        damage_ax.grid(True)


def padded_limits(values: list[np.ndarray]) -> tuple[float, float]:
    lo = min(float(np.nanmin(v)) for v in values)
    hi = max(float(np.nanmax(v)) for v in values)
    if np.isclose(lo, hi):
        pad = max(abs(lo) * 0.1, 1e-12)
    else:
        pad = 0.05 * (hi - lo)
    return lo - pad, hi + pad


def save_movie(files: list[Path], args: argparse.Namespace) -> None:
    frames = [(conf_index(path), load_conf(path)) for path in files]
    ncols = 1 + int(args.strain_zoom) + int(args.displacement) + int(args.damage)
    fig, axes = make_axes(ncols)
    strain_ax, zoom_ax, disp_ax, damage_ax = selected_axes(
        axes, args.displacement, args.damage, args.strain_zoom
    )
    style_axes(strain_ax, zoom_ax, disp_ax, damage_ax)

    x_values = [data[:, 0] for _, data in frames]
    strain_values = [data[:, 2] for _, data in frames]
    strain_ax.set_xlim(*padded_limits(x_values))
    strain_ax.set_ylim(*padded_limits(strain_values))
    (strain_line,) = strain_ax.plot([], [], lw=2)

    zoom_line = None
    if zoom_ax is not None:
        zoom_mask = [
            (data[:, 0] >= args.zoom_x_min) & (data[:, 0] <= args.zoom_x_max)
            for _, data in frames
        ]
        zoom_values = [
            data[:, 2][mask] for (_, data), mask in zip(frames, zoom_mask) if np.any(mask)
        ]
        zoom_ax.set_xlim(args.zoom_x_min, args.zoom_x_max)
        zoom_ax.set_ylim(*padded_limits(zoom_values))
        (zoom_line,) = zoom_ax.plot([], [], lw=2)

    disp_line = None
    if disp_ax is not None:
        disp_values = [data[:, 1] for _, data in frames]
        disp_ax.set_xlim(*padded_limits(x_values))
        disp_ax.set_ylim(*padded_limits(disp_values))
        (disp_line,) = disp_ax.plot([], [], lw=2)

    damage_line = None
    if damage_ax is not None:
        damage_ax.set_xlim(*padded_limits(x_values))
        (damage_line,) = damage_ax.plot([], [], lw=2)

    fig.tight_layout()
    writer = FFMpegWriter(fps=args.fps, metadata={"title": "Peridynamics profiles"})
    with writer.saving(fig, str(args.movie), dpi=200):
        for step, data in frames:
            x = data[:, 0]
            strain_line.set_data(x, data[:, 2])
            if zoom_line is not None:
                mask = (x >= args.zoom_x_min) & (x <= args.zoom_x_max)
                zoom_line.set_data(x[mask], data[:, 2][mask])
                if args.dynamic_zoom_y:
                    zoom_ax.set_ylim(*padded_limits([data[:, 2][mask]]))
            if disp_line is not None:
                disp_line.set_data(x, data[:, 1])
            if damage_line is not None:
                damage_line.set_data(x, data[:, 3])
            fig.suptitle(f"Increment {step}", y=0.995)
            writer.grab_frame()
    plt.close(fig)


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
        "--strain-zoom",
        action="store_true",
        help="Also plot a zoomed strain panel.",
    )
    parser.add_argument(
        "--zoom-x-min",
        default=0.05,
        type=float,
        help="Minimum x-position for --strain-zoom.",
    )
    parser.add_argument(
        "--zoom-x-max",
        default=0.95,
        type=float,
        help="Maximum x-position for --strain-zoom.",
    )
    parser.add_argument(
        "--dynamic-zoom-y",
        action="store_true",
        help="Rescale the strain zoom y-axis on each movie frame.",
    )
    parser.add_argument(
        "--save",
        type=Path,
        help="Save figure to this file instead of only showing it.",
    )
    parser.add_argument(
        "--movie",
        type=Path,
        help="Save an MP4 movie from the selected conf*.dat files.",
    )
    parser.add_argument(
        "--fps",
        default=12,
        type=int,
        help="Frames per second for --movie output.",
    )
    parser.add_argument(
        "--no-legend",
        action="store_true",
        help="Do not draw legends on static plots.",
    )
    args = parser.parse_args()

    files = sorted(args.result_dir.glob("conf*.dat"), key=conf_index)
    if args.steps:
        wanted = set(args.steps)
        files = [path for path in files if conf_index(path) in wanted]

    if not files:
        raise SystemExit(f"No conf*.dat files found in {args.result_dir}")

    if args.movie:
        save_movie(files, args)
        return

    ncols = 1 + int(args.strain_zoom) + int(args.displacement) + int(args.damage)
    fig, axes = make_axes(ncols)
    strain_ax, zoom_ax, disp_ax, damage_ax = selected_axes(
        axes, args.displacement, args.damage, args.strain_zoom
    )

    for path in files:
        data = load_conf(path)
        x = data[:, 0]
        u = data[:, 1]
        strain = data[:, 2]
        damage = data[:, 3]
        label = f"conf {conf_index(path)}"

        strain_ax.plot(x, strain, label=label)
        if zoom_ax is not None:
            mask = (x >= args.zoom_x_min) & (x <= args.zoom_x_max)
            zoom_ax.plot(x[mask], strain[mask], label=label)
        if disp_ax is not None:
            disp_ax.plot(x, u, label=label)
        if damage_ax is not None:
            damage_ax.plot(x, damage, label=label)

    if zoom_ax is not None:
        zoom_ax.set_xlim(args.zoom_x_min, args.zoom_x_max)
    style_axes(strain_ax, zoom_ax, disp_ax, damage_ax)
    if not args.no_legend:
        strain_ax.legend()
        if zoom_ax is not None:
            zoom_ax.legend()
        if disp_ax is not None:
            disp_ax.legend()
        if damage_ax is not None:
            damage_ax.legend()

    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=200)
    plt.show()


if __name__ == "__main__":
    main()
