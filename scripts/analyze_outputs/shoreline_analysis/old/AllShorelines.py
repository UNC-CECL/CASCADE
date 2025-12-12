import os
import argparse
import numpy as np
import matplotlib.pyplot as plt


def load_cascade(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    return data["cascade"][0]


def get_x_s_TS(b3d):
    """CASCADE sometimes stores shoreline as x_s_TS or _x_s_TS."""
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    elif hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    else:
        raise AttributeError("No shoreline time series found on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True, relative=False):
    """
    Build shoreline[t, domain].
    """
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)

    # Time dimension from domain 0
    nt = len(get_x_s_TS(b3d_list[0]))

    shoreline = np.zeros((nt, ndom))

    for j in range(ndom):
        xs = get_x_s_TS(b3d_list[j])  # dam

        shoreline[:, j] = xs

    return shoreline


def plot_shoreline_snapshots_by_domain(shoreline, outdir, start_year=0):
    """
    For each time step (year), plot shoreline vs domain number.
    """

    os.makedirs(outdir, exist_ok=True)

    nt, ndom = shoreline.shape
    domain_numbers = np.arange(ndom)

    for t in range(nt):
        year = start_year + t

        plt.figure(figsize=(10, 4))
        plt.plot(domain_numbers, shoreline[t, :], marker="o", linewidth=1.2)

        plt.xlabel("Domain number")
        plt.ylabel("Shoreline position (m)")
        plt.title(f"Shoreline Position Across Domains – Year {year}")
        plt.grid(alpha=0.3)
        plt.tight_layout()

        fname = os.path.join(outdir, f"shoreline_snapshot_year_{year}.png")
        plt.savefig(fname, dpi=200)
        plt.close()

        print(f"Saved: {fname}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("npz_file", help="Path to CASCADE .npz file")
    parser.add_argument(
        "--outdir",
        default="shoreline_snapshots",
        help="Output directory"
    )
    parser.add_argument(
        "--start_year",
        type=int,
        default=0,
        help="Label for starting year"
    )
    parser.add_argument(
        "--relative",
        action="store_true",
        help="Plot shoreline change relative to t=0 instead of absolute position"
    )
    args = parser.parse_args()

    cascade = load_cascade(args.npz_file)

    shoreline = build_shoreline_matrix(
        cascade,
        to_meters=True,
        relative=args.relative,
    )

    plot_shoreline_snapshots_by_domain(
        shoreline,
        outdir=args.outdir,
        start_year=args.start_year,
    )


if __name__ == "__main__":
    main()
