import os

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker as mticker
import numpy as np

from ..vtk import VTKVertexReader
from ..fibergrowth.evaluate import Fiber

from ._utils import load_iteration_results, load_fibers


def get_filepaths_to_plot(df, sample_size):
    yaml_paths = df["filepath"].iloc[:sample_size]
    vtk_paths = [os.path.splitext(path)[0] + ".vtu" for path in yaml_paths]
    return yaml_paths, vtk_paths


def plot_mean_displacement(
    vtk_files,
    yaml_cfg_files,
    outfile,
    displacement_dataset="Displacement Field_0_",
    stress_dataset="vonmises_",
    tri_subdiv=3,
):
    """Plot the mean of a set of simulation outputs.

    This computes the mean displacement and stress fields and plots a visualization
    into the specified output file.
    """
    # Compute mean fields
    displacement_mean = []
    stress_mean = []
    for vtk_file in vtk_files:
        vtk_reader = VTKVertexReader(vtk_file)
        displacement_mean.append(vtk_reader[displacement_dataset])
        stress_mean.append(vtk_reader[stress_dataset])
    displacement_mean = np.mean(displacement_mean, axis=0)
    stress_mean = np.mean(stress_mean, axis=0)

    # Take points and connectivity from last file
    # NOTE: This assumes that this information is the same for all files
    #       (which is not checked here)!
    points = vtk_reader.points
    points_displ = points + displacement_mean
    connect = vtk_reader.connectivity

    # Triangulate, interpolate, and refine
    tri = mtri.Triangulation(points_displ[..., 0], points_displ[..., 1], connect)
    stress_interp = mtri.LinearTriInterpolator(tri, stress_mean)
    refiner = mtri.UniformTriRefiner(tri)
    tri_ref, stress_ref = refiner.refine_field(
        stress_mean, subdiv=tri_subdiv, triinterpolator=stress_interp
    )

    # Single list of all fibers
    fibers = sum([load_fibers(file) for file in yaml_cfg_files], start=[])

    # Create figure
    fig, ax = plt.subplots(1, 1, constrained_layout=True, dpi=300)
    ax.grid(True, zorder=-10)
    locator = mticker.MaxNLocator()

    # Plot stress and grid
    cb_data = ax.tricontourf(
        tri_ref, stress_ref, locator=locator, zorder=2, cmap="YlOrRd"
    )
    # ax.triplot(tri, "-", linewidth=0.5, color="k", markersize=1.0, zorder=3)

    # Plot fibers
    tri = mtri.Triangulation(points[..., 0], points[..., 1], connect)
    displ_interp_x = mtri.LinearTriInterpolator(tri, displacement_mean[..., 0])
    displ_interp_y = mtri.LinearTriInterpolator(tri, displacement_mean[..., 1])
    for fiber in fibers:
        x_vals = np.linspace(fiber["start"][0], fiber["end"][0], 100)
        y_vals = np.linspace(fiber["start"][1], fiber["end"][1], 100)
        x_vals_displ = x_vals + displ_interp_x(x_vals, y_vals)
        y_vals_displ = y_vals + displ_interp_y(x_vals, y_vals)
        ax.plot(x_vals_displ, y_vals_displ, color="k", zorder=5, alpha=0.1)

    # Colorbar
    cb = plt.colorbar(cb_data, ax=ax, orientation="vertical")
    cb.set_label(r"Von Mises Stress $[\mathrm{Pa}]$")

    # Figure aesthetics
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
        description="DUNE STRUCTURES Population Mean Plot CLI"
    )
    parser.add_argument(
        "root_dir",
        type=os.path.realpath,
        help="Root directory of the optimization algorithm output",
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
    # parser.add_argument("--plot-fittest", "-f", action="store_true")
    # parser.add_argument("--plot-lowest-stress", "-s", action="store_true")
    # parser.add_argument("--plot-lowest-fiber", "-l", action="store_true")
    args = parser.parse_args()

    # Load data
    print("Loading CSV data...")
    results = load_iteration_results(args.root_dir)

    # Plot population means
    print("Plotting means for {} iterations...".format(len(results)))
    for i, df in enumerate(results):
        for data_plot, outname in [
            (df.sort_values(["fitness", "stress"], ascending=[False, True]), "fittest"),
            (df.sort_values("stress", ascending=True), "stress-low"),
            (df.sort_values("length", ascending=True), "fiber-low"),
        ]:
            yaml_paths, vtk_paths = get_filepaths_to_plot(data_plot, args.sample_size)
            plot_mean_displacement(
                vtk_paths,
                yaml_paths,
                os.path.join(
                    args.root_dir,
                    "{}-{:03d}.{}".format(outname, i, args.extension),
                ),
            )
