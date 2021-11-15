import os

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker as mticker
import numpy as np

from ruamel.yaml import YAML

from ..vtk import VTKVertexReader
from ..fibergrowth.evaluate import Fiber


def load_yaml_config(filename):
    with open(filename, "r") as file:
        yaml = YAML(typ="safe")
        return yaml.load(file)


def get_block_by_name_recursive(yaml_node, blockname):
    for block in yaml_node:
        if block["_blockname"] == blockname:
            return block
        elif "blocks" in block:
            return get_block_by_name_recursive(block["blocks"], blockname)
    return None


def plot_mean_displacement(
    vtk_files, yaml_cfg_files, outfile, displacement_dataset, stress_dataset, tri_subdiv
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

    # Collect fiber data
    def fetch_fibers(yaml_file_path):
        cfg = load_yaml_config(yaml_file_path)
        return [
            fiber
            for fiber in get_block_by_name_recursive(
                cfg["solver"]["blocks"], "linearsolver_0"
            )["operator"]["reinforced_operator"]["fibres"]
        ]

    # Single list of all fibers
    fibers = sum([fetch_fibers(file) for file in yaml_cfg_files], start=[])

    # Create figure
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    ax.grid(True, zorder=-10)
    locator = mticker.MaxNLocator()

    # Plot stress and grid
    cb_data = ax.tricontourf(tri_ref, stress_ref, locator=locator, zorder=2)
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
    cb.set_label(r"Von Mises Stress $\sigma_v / \mathrm{Pa}$")

    # Figure aesthetics
    ax.set_aspect("equal", anchor="C", adjustable="datalim")
    ax.set_xlabel(r"Extension $x / \mathrm{m}$")
    ax.set_ylabel(r"Extension $y / \mathrm{m}$")

    # Save file
    if not os.path.splitext(outfile)[1]:
        outfile += ".pdf"
    fig.savefig(outfile)
    plt.close(fig)
