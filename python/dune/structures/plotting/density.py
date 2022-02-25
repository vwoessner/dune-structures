import os
import glob

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker as mticker
import matplotlib.collections as mcol
import numpy as np

from ..vtk import VTKVertexReader
from ..fibergrowth.evaluate import Fiber

from ._utils import load_iteration_results, load_fibers, get_filepaths_to_plot
from .displacement import load_data


def plot_fiber_density(
    vtk_files,
    yaml_cfg_files,
    outfile,
    displacement_dataset="Displacement Field_0_",
    prec=100,
):
    # Compute mean displacement
    displacement_mean = np.mean(
        [VTKVertexReader(vtk_file)[displacement_dataset] for vtk_file in vtk_files],
        axis=0,
    )

    # Collect all fibers
    fibers = sum([load_fibers(file) for file in yaml_cfg_files], start=[])

    # Compute displaced geometry
    vtk_reader = VTKVertexReader(vtk_files[0])
    points = vtk_reader.points
    connect = vtk_reader.connectivity
    points_displ = points + displacement_mean
    points_displ = points_displ[..., :2]

    # Compute reference geometry
    tri = mtri.Triangulation(points[..., 0], points[..., 1], connect)
    trifinder = tri.get_trifinder()

    # Compute fiber density for every triangle
    density = np.zeros((len(connect)), dtype=np.int64)
    for fiber in fibers:
        # Store all triangles this fiber crosses
        x_vals = np.linspace(fiber["start"][0], fiber["end"][0], prec)
        y_vals = np.linspace(fiber["start"][1], fiber["end"][1], prec)
        triangles = set(trifinder(x_vals, y_vals))
        triangles.discard(-1)
        triangles = list(triangles)
        density[triangles] = density[triangles] + 1
    density = np.true_divide(density, np.max(density), dtype=np.float64)

    # Create Polygon collection
    polys = mcol.PolyCollection(
        [points_displ[c] for c in connect], closed=False, cmap="viridis",
        zorder=3
    )
    polys.set_array(density)

    # Create figure
    fig, ax = plt.subplots(1, 1, constrained_layout=True, dpi=300)
    ax.grid(True, zorder=-10)

    # Plot polygons
    ax.add_collection(polys)

    # Colorbar
    cb = fig.colorbar(polys, ax=ax, orientation="vertical")
    cb.set_label("Relative Stress Fiber Density [–]")

    # Aesthetics
    ax.autoscale_view()
    ax.set_aspect("equal", anchor="C", adjustable="datalim")
    ax.set_xlabel(r"Extension $[\mathrm{\mu m}]$")
    ax.set_ylabel(r"Extension $[\mathrm{\mu m}]$")

    # Save file
    if not os.path.splitext(outfile)[1]:
        outfile += ".pdf"
    fig.savefig(outfile)
    plt.close(fig)


def entrypoint():
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(
        description="DUNE STRUCTURES Stress Fiber Density CLI"
    )
    parser.add_argument(
        "root_dir",
        type=os.path.realpath,
        help="Root directory of the optimization algorithm output",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=os.path.relpath,
        help="Output directory relative to root_dir",
        default="fiber-densities",
    )
    parser.add_argument(
        "--extension",
        "-e",
        type=str,
        help="Output file extension (without dot)",
        default="pdf",
    )
    parser.add_argument(
        "--sample-size",
        "-n",
        type=int,
        help="Number of samples taken to compute the mean",
        default=20,
    )
    parser.add_argument(
        "--precision",
        "-p",
        type=int,
        help="Fiber sampling precision",
        default=200,
    )

    args = parser.parse_args()

    print("Loading CSV data...")
    results = load_iteration_results(args.root_dir)

    # Plot densities
    outdir = os.path.join(os.path.abspath(args.root_dir), args.output_dir)
    os.makedirs(outdir)

    print(
        "Plotting fiber densities for {} iterations into {}...".format(
            len(results), outdir
        )
    )
    for i, df in enumerate(results):
        for data_plot, outname in [
            (df.sort_values(["fitness", "stress"], ascending=[False, True]), "fittest"),
            (df.sort_values("stress", ascending=True), "stress-low"),
            (df.sort_values("length", ascending=True), "fiber-low"),
        ]:
            yaml_paths, vtk_paths = get_filepaths_to_plot(data_plot, args.sample_size)
            plot_fiber_density(
                vtk_paths,
                yaml_paths,
                os.path.join(
                    outdir,
                    "{}-{:03d}.{}".format(outname, i, args.extension),
                ),
                prec=args.precision
            )
