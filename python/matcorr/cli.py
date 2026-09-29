"""
Command-line interface for matcorr (Materials Correlation Analysis).

Provides headless terminal-based structural analysis for liquid, amorphous,
and nanostructured materials simulations.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import correlation
from matcorr import __version__


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments for matcorr CLI."""
    parser = argparse.ArgumentParser(
        prog="matcorr",
        description="matcorr — High-performance structural analysis for atomistic simulations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        nargs="?",
        help="Input trajectory or structure file (XYZ, CAR, CIF, XDATCAR, GRO, etc.)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_base",
        default=None,
        help="Output base path for generated files (default: input file stem)",
    )
    parser.add_argument(
        "-m",
        "--material",
        choices=["amorphous", "liquid", "crystalline"],
        default="amorphous",
        help="Material preset: automatically configures optimal binning and smoothing",
    )
    parser.add_argument(
        "--r-max",
        type=float,
        default=20.0,
        help="Maximum radius for pair/radial distribution function (Angstroms)",
    )
    parser.add_argument(
        "--r-bin",
        type=float,
        default=None,
        help="Bin width for radial distribution function (default: preset-dependent)",
    )
    parser.add_argument(
        "--r-int-max",
        type=float,
        default=10.0,
        help="Maximum radius for coordination number integration",
    )
    parser.add_argument(
        "--q-max",
        type=float,
        default=20.0,
        help="Maximum scattering vector magnitude Q for S(Q) (1/Angstroms)",
    )
    parser.add_argument(
        "--q-bin",
        type=float,
        default=None,
        help="Bin width for structure factor S(Q) (default: preset-dependent)",
    )
    parser.add_argument(
        "--angle-bin",
        type=float,
        default=None,
        help="Bin width for planar angle distribution in degrees (default: preset-dependent)",
    )
    parser.add_argument(
        "--dihedral-bin",
        type=float,
        default=None,
        help="Bin width for dihedral angle distribution in degrees",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=3.0,
        help="Bond cutoff radius for neighbor searching (Angstroms)",
    )
    parser.add_argument(
        "--smoothing-sigma",
        type=float,
        default=None,
        help="Bandwidth sigma for post-processing kernel smoothing",
    )
    parser.add_argument(
        "--smoothing-kernel",
        choices=["gaussian", "bump", "triweight", "epanechnikov", "cosine", "biweight"],
        default="gaussian",
        help="Kernel smoothing function",
    )
    parser.add_argument(
        "--no-smoothing",
        action="store_true",
        help="Disable post-processing kernel smoothing",
    )
    parser.add_argument(
        "--csv",
        dest="csv",
        action="store_true",
        default=True,
        help="Export results to CSV files (enabled by default)",
    )
    parser.add_argument(
        "--no-csv",
        dest="csv",
        action="store_false",
        help="Disable CSV file export",
    )
    parser.add_argument(
        "--min-frame",
        type=int,
        default=1,
        help="1-based starting trajectory frame index",
    )
    parser.add_argument(
        "--max-frame",
        type=int,
        default=-1,
        help="1-based ending trajectory frame index (-1 for all frames)",
    )
    parser.add_argument(
        "--time-step",
        type=float,
        default=1.0,
        help="Simulation time step in femtoseconds",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress console progress and telemetry",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"matcorr {__version__} (Correlation analysis suite)",
    )

    return parser.parse_args(argv)


def apply_material_defaults(args: argparse.Namespace) -> None:
    """Apply default physical binning parameters based on material type."""
    if args.material == "crystalline":
        if args.r_bin is None:
            args.r_bin = 0.002
        if args.q_bin is None:
            args.q_bin = 0.002
        if args.angle_bin is None:
            args.angle_bin = 0.1
        if args.smoothing_sigma is None:
            args.smoothing_sigma = 0.01
    elif args.material == "liquid":
        if args.r_bin is None:
            args.r_bin = 0.05
        if args.q_bin is None:
            args.q_bin = 0.05
        if args.angle_bin is None:
            args.angle_bin = 0.5
        if args.smoothing_sigma is None:
            args.smoothing_sigma = 0.15
    else:  # amorphous
        if args.r_bin is None:
            args.r_bin = 0.02
        if args.q_bin is None:
            args.q_bin = 0.02
        if args.angle_bin is None:
            args.angle_bin = 1.0
        if args.smoothing_sigma is None:
            args.smoothing_sigma = 0.1

    if args.dihedral_bin is None:
        args.dihedral_bin = args.angle_bin


def main(argv: list[str] | None = None) -> int:
    """Main CLI entrypoint for matcorr."""
    args = parse_args(argv)

    if not args.input_file:
        print("matcorr: error: the following arguments are required: input_file", file=sys.stderr)
        print("Run 'matcorr --help' for usage details.", file=sys.stderr)
        return 1

    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"matcorr: error: file not found: {input_path}", file=sys.stderr)
        return 1

    apply_material_defaults(args)

    output_base = args.output_base or str(input_path.with_suffix(""))
    output_dir = Path(output_base).parent
    if output_dir and not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    if not args.quiet:
        print("─" * 60)
        print(f"  matcorr {__version__} — Materials Correlation Engine")
        print("─" * 60)
        print(f"  Input File  : {input_path}")
        print(f"  Output Base : {output_base}")
        print(f"  Material    : {args.material}")
        print(f"  r-range     : [0.0, {args.r_max:.2f}] Å (bin: {args.r_bin:.4f} Å)")
        print(f"  q-range     : [0.0, {args.q_max:.2f}] Å⁻¹ (bin: {args.q_bin:.4f} Å⁻¹)")
        print(f"  Angle bin   : {args.angle_bin:.2f}°")
        print(f"  Smoothing   : {'OFF' if args.no_smoothing else f'{args.smoothing_kernel} (σ={args.smoothing_sigma:.4f})'}")
        print("─" * 60)

    start_time = time.perf_counter()

    try:
        data = correlation.read(str(input_path))
    except Exception as exc:
        print(f"matcorr: failed to read input file: {exc}", file=sys.stderr)
        return 1

    # Single Cell vs Trajectory handling
    if isinstance(data, correlation.Cell):
        cell = data
        num_atoms = cell.atom_count() if callable(cell.atom_count) else cell.atom_count
        if not args.quiet:
            print(f"  Single frame loaded: {num_atoms} atoms")

        df = correlation.DistributionFunctions(cell, cutoff=args.cutoff)
        df.calculate_rdf(r_max=args.r_max, bin_width=args.r_bin)

        try:
            df.calculate_pad(bin_width=args.angle_bin)
        except Exception:
            pass

        if not args.no_smoothing:
            try:
                kernel_map = {
                    "gaussian": correlation.KernelType.Gaussian,
                    "bump": correlation.KernelType.Bump,
                    "triweight": correlation.KernelType.Triweight,
                }
                k_type = kernel_map.get(args.smoothing_kernel, correlation.KernelType.Gaussian)
                df.smooth_all(sigma=args.smoothing_sigma, kernel=k_type)
            except Exception:
                pass

        if args.csv:
            correlation.write_csv(output_base, df, write_smoothed=not args.no_smoothing)
            if not args.quiet:
                print(f"  Results saved to CSV: {output_base}_*.csv")

    elif isinstance(data, correlation.Trajectory):
        traj = data
        n_frames = traj.num_frames()
        if not args.quiet:
            print(f"  Trajectory loaded: {n_frames} frames")

        # Process frames within requested bounds
        start_f = max(0, args.min_frame - 1)
        end_f = n_frames if args.max_frame < 0 else min(n_frames, args.max_frame)

        if start_f >= end_f:
            print(f"matcorr: error: invalid frame range [{args.min_frame}, {args.max_frame}]", file=sys.stderr)
            return 1

        first_cell = traj[start_f]
        accum_df = correlation.DistributionFunctions(first_cell, cutoff=args.cutoff)
        accum_df.calculate_rdf(r_max=args.r_max, bin_width=args.r_bin)

        for frame_idx in range(start_f + 1, end_f):
            c = traj[frame_idx]
            frame_df = correlation.DistributionFunctions(c, cutoff=args.cutoff)
            frame_df.calculate_rdf(r_max=args.r_max, bin_width=args.r_bin)
            accum_df.add(frame_df)

        if not args.no_smoothing:
            try:
                kernel_map = {
                    "gaussian": correlation.KernelType.Gaussian,
                    "bump": correlation.KernelType.Bump,
                    "triweight": correlation.KernelType.Triweight,
                }
                k_type = kernel_map.get(args.smoothing_kernel, correlation.KernelType.Gaussian)
                accum_df.smooth_all(sigma=args.smoothing_sigma, kernel=k_type)
            except Exception:
                pass

        if args.csv:
            correlation.write_csv(output_base, accum_df, write_smoothed=not args.no_smoothing)
            if not args.quiet:
                print(f"  Averaged results saved to CSV: {output_base}_*.csv")

    else:
        print("matcorr: error: unrecognized data object returned by reader", file=sys.stderr)
        return 1

    elapsed = time.perf_counter() - start_time
    if not args.quiet:
        print("─" * 60)
        print(f"  Analysis completed successfully in {elapsed:.3f} s")
        print("─" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
