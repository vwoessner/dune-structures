import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.collections as mcollections
import matplotlib
import numpy as np
import meshio


def plot_displacement(
    ax,
    vmin,
    vmax,
    points,
    connectivity,
    displacement,
    stress,
    fibers,
    plot_grid=True,
    fiber_alpha=1.0,
    plot_fiber=True,
):
    # Triangulate, interpolate, and refine
    points_displ = points + displacement
    tri = mtri.Triangulation(points_displ[..., 0], points_displ[..., 1], connectivity)
    # stress_interp = mtri.LinearTriInterpolator(tri, stress)
    # refiner = mtri.UniformTriRefiner(tri)
    # tri_ref, stress_ref = refiner.refine_field(
    #     stress, subdiv=tri_subdiv, triinterpolator=stress_interp
    # )

    # Create figure

    # locator = mticker.MaxNLocator()

    cb_data = ax.tripcolor(
        tri,
        stress,
        #vmin=vmin,
        #vmax=vmax,
        zorder=2,
        cmap="YlOrRd",
        rasterized=True,
        shading="flat",
        # cmap=tol_cmap("iridescent"),
        norm=mcolors.LogNorm(vmin=vmin, vmax=vmax),
    )
    if plot_grid:
        ax.triplot(
            tri, "-", linewidth=0.5, color="k", markersize=1.0, zorder=3, alpha=0.8
        )

    ax.grid(True, zorder=-10)

    if plot_fiber==True:
        tri = mtri.Triangulation(points[..., 0], points[..., 1], connectivity)
        displ_interp_x = mtri.LinearTriInterpolator(tri, displacement[..., 0])
        displ_interp_y = mtri.LinearTriInterpolator(tri, displacement[..., 1])
        for fiber in fibers:
            x_vals = np.linspace(fiber["start"][0], fiber["end"][0], 100)
            y_vals = np.linspace(fiber["start"][1], fiber["end"][1], 100)
            x_vals_displ = x_vals + displ_interp_x(x_vals, y_vals)
            y_vals_displ = y_vals + displ_interp_y(x_vals, y_vals)
            ax.plot(x_vals_displ, y_vals_displ, color="k", zorder=5, alpha=fiber_alpha, linewidth=1.0)

    return cb_data
    
    
def load_mesh(meshfile):
    mesh = meshio.read(meshfile)
    # NOTE: Meshio v3.3.1
    # NOTE: Assuming triangles
    return {"points": mesh.points, "connectivity": mesh.cells["triangle"]}
    
    
def plot_mesh(
    ax,
    meshfile,
    fibers,
    plot_grid=False,
    fiber_alpha=1.0,
    plot_fiber=True
    ):
    mesh = load_mesh(meshfile)
    points = mesh["points"]
    connectivity = mesh["connectivity"]
    
    points_displ = points
    x = points[...,0]
    y = points[...,1]

    # Tripcolor plot
    tri = mtri.Triangulation(points_displ[..., 0], points_displ[..., 1], connectivity)
    # Use binary map, such that Matplotlib-Grid is not visible
    cb_data = ax.tripcolor(tri, facecolors=np.zeros((len(connectivity))), cmap="binary", linewidth=0.3, edgecolor="k", zorder=3, alpha=1.0)

    ax.grid(True, zorder=-10)
    
    for fiber in fibers:
        if fiber["end"][1] > np.max(points[...,1]):
            fiber["end"][1] = np.max(points[...,1])
        x_vals = np.linspace(fiber["start"][0], fiber["end"][0], 100, endpoint=True)
        y_vals = np.linspace(fiber["start"][1], fiber["end"][1], 100, endpoint=True)
        x_vals_displ = x_vals
        y_vals_displ = y_vals
        if plot_fiber==True:
            ax.plot(x_vals_displ, y_vals_displ, color="k", zorder=5, alpha=fiber_alpha, linewidth=1.0)
            
    if plot_grid:
        ax.triplot(tri, "-", linewidth=0.5, color="k", markersize=1.0, zorder=3, alpha=0.8)

    return cb_data
