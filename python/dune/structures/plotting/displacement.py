import os
import glob

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker as mticker
import numpy as np

from ..vtk import VTKVertexReader
from ..fibergrowth.evaluate import Fiber

from ._utils import load_iteration_results, load_fibers, get_filepaths_to_plot

# plt.style.use('dark_background')

def load_data(
    vtk_file,
    yaml_input_file,
    displacement_dataset="Displacement Field_0_",
    stress_dataset="vonmises_",
):
    fibers = load_fibers(yaml_input_file)
    vtk_reader = VTKVertexReader(vtk_file)
    return {
        "fibers": fibers,
        "points": vtk_reader.points,
        "connectivity": vtk_reader.connectivity,
        "stress": vtk_reader[stress_dataset],
        "displacement": vtk_reader[displacement_dataset],
    }


def plot_displacement(
    points,
    connectivity,
    displacement,
    stress,
    fibers,
    outfile,
    tri_subdiv=3,
    plot_grid=True,
    fiber_alpha=1.0,
):
    # Triangulate, interpolate, and refine
    points_displ = points + displacement
    tri = mtri.Triangulation(points_displ[..., 0], points_displ[..., 1], connectivity)
    stress_interp = mtri.LinearTriInterpolator(tri, stress)
    refiner = mtri.UniformTriRefiner(tri)
    tri_ref, stress_ref = refiner.refine_field(
        stress, subdiv=tri_subdiv, triinterpolator=stress_interp
    )

    # Create figure
    fig, ax = plt.subplots(1, 1, constrained_layout=True, dpi=300)
    ax.grid(True, zorder=-10)
    locator = mticker.MaxNLocator()

    cb_data = ax.tricontourf(
        tri_ref, stress_ref, locator=locator, zorder=2, cmap="YlOrRd"
    )
    if plot_grid:
        ax.triplot(
            tri, "-", linewidth=0.5, color="k", markersize=1.0, zorder=3, alpha=0.8
        )

    tri = mtri.Triangulation(points[..., 0], points[..., 1], connectivity)
    displ_interp_x = mtri.LinearTriInterpolator(tri, displacement[..., 0])
    displ_interp_y = mtri.LinearTriInterpolator(tri, displacement[..., 1])
    for fiber in fibers:
        x_vals = np.linspace(fiber["start"][0], fiber["end"][0], 100)
        y_vals = np.linspace(fiber["start"][1], fiber["end"][1], 100)
        x_vals_displ = x_vals + displ_interp_x(x_vals, y_vals)
        y_vals_displ = y_vals + displ_interp_y(x_vals, y_vals)
        ax.plot(x_vals_displ, y_vals_displ, color="k", zorder=5, alpha=fiber_alpha)

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
    fibers = sum([load_fibers(file) for file in yaml_cfg_files], start=[])

    # Take points and connectivity from last file
    # NOTE: This assumes that this information is the same for all files
    #       (which is not checked here)!
    plot_displacement(
        vtk_reader.points,
        vtk_reader.connectivity,
        displacement_mean,
        stress_mean,
        fibers,
        outfile,
        tri_subdiv,
        plot_grid=False,
        fiber_alpha=0.1,
    )


def entrypoint_single():
    import argparse

    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Plot displacement and stress of a 2D domain from a VTK file"
    )
    parser.add_argument("input_file", help="Input VTK file", type=os.path.realpath)
    parser.add_argument(
        "--dataset-displacement",
        help="Name of the dataset containing displacement data",
        type=str,
        default="Displacement Field_0_",
    )
    parser.add_argument(
        "--dataset-stress",
        help="Name of the dataset containing stress data",
        type=str,
        default="vonmises_",
    )
    parser.add_argument(
        "--output-file",
        "-o",
        help="Path of the plot output file",
        type=os.path.realpath,
        default="plot.pdf",
    )
    parser.add_argument(
        "--yaml-input-file",
        "-i",
        help="Input file of simulation containing fiber location data",
        type=os.path.realpath,
    )
    parser.add_argument(
        "--precision",
        help="Output plotting precision per standard length",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--contour-ref", help="Refinement for the contour plot", type=int, default=3
    )
    parser.add_argument(
        "--glob",
        help="Input file contains wildcard. Infer input and output file names from input name",
        action="store_true",
    )
    args = parser.parse_args()

    if not args.glob:
        plot_displacement(
            outfile=args.output_file,
            tri_subdiv=args.contour_ref,
            **load_data(
                args.input_file,
                args.yaml_input_file,
                args.dataset_displacement,
                args.dataset_stress,
            )
        )
    # Special handling for multiple input files
    else:
        for file in sorted(glob.glob(args.input_file)):
            try:
                filepath, _ = os.path.splitext(file)
                # Manipulate stored arguments
                args.input_file = file
                args.yaml_input_file = filepath + ".yml"
                args.output_file = filepath + ".png"
                plot_displacement(
                    outfile=args.output_file,
                    tri_subdiv=args.contour_ref,
                    **load_data(
                        args.input_file,
                        args.yaml_input_file,
                        args.dataset_displacement,
                        args.dataset_stress,
                    )
                )
            except RuntimeError as e:
                import warnings

                warnings.warn("Unable to visualize file: {}".format(file))


def entrypoint_mean():
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
