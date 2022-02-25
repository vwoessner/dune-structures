import os

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import matplotlib.lines as mlines
import pandas as pd
import numpy as np

from ._utils import load_iteration_results


def _limits_dict():
    return {"stress": {}, "volume": {}, "length": {}}


def determine_plot_limits(iteration_results):
    volume_min = np.min([df["length"].min() for df in iteration_results])
    volume_max = np.max([df["length"].max() for df in iteration_results])
    stress_min = np.min([df["stress"].min() for df in iteration_results])
    stress_max = np.median([df["stress"].max() for df in iteration_results])

    def fitness(df):
        return df["fitness"].loc[df["selected"]]

    fitness_min = np.min([fitness(df).min() for df in iteration_results])
    fitness_max = np.max([fitness(df).max() for df in iteration_results])

    return {
        "stress": {"min": stress_min, "max": stress_max},
        "volume": {"min": volume_min, "max": volume_max},
        "fitness": {"min": fitness_min, "max": fitness_max},
    }


def plot_population(df_iteration, iteration, outfile, limits=None):
    # Extract the data
    select = df_iteration["selected"]
    volume = df_iteration["length"]
    stress = df_iteration["stress"]
    this_it = df_iteration["iteration"] == iteration

    # Create the figure
    fig, ax = plt.subplots(1, 1, constrained_layout=True, dpi=300)

    # If limits are given, draw the points transparently and use autoscale on them
    if limits is not None:
        ax.plot(
            [limits["volume"]["min"], limits["volume"]["max"]],
            [limits["stress"]["min"], limits["stress"]["max"]],
            "o",
            alpha=0.0,
        )
        ax.autoscale_view()
        ax.autoscale(False)

    # Plot removed
    red = mcm.get_cmap("Reds")(0.7)
    ax.plot(
        volume.loc[~select],
        stress.loc[~select],
        "o",
        label="Removed",
        markeredgecolor=red,
        fillstyle="none",
    )

    # Plot selected
    fitness = df_iteration["fitness"].loc[select]
    # Use calculated or current limits
    lim_fitness = limits["fitness"] if limits else {}
    norm = mcolors.Normalize(
        lim_fitness.get("min", fitness.min())
        - 0.2
        * (
            lim_fitness.get("max", fitness.max())
            - lim_fitness.get("min", fitness.min())
        ),
        lim_fitness.get("max", fitness.max()),
    )
    cb_data = ax.scatter(
        volume[select],
        stress[select],
        c=fitness,
        cmap="Blues",
        label="Selected",
        norm=norm,
    )

    # Plot created
    ax.plot(
        volume.loc[this_it],
        stress.loc[this_it],
        "x",
        markersize=5,
        label="Created",
        color="grey",
    )

    # Colorbar
    cb = fig.colorbar(cb_data, ax=ax)
    cb.set_label("Relative Fitness")
    cb.ax.set_ylim(bottom=fitness.min())  # Do this in case limits are not set

    # Plot mean values
    # ax.axhline(np.mean(stresses[~select]), color="r")
    # ax.axhline(np.mean(stresses[select]), color="b")

    # Aesthetics
    ax.set_title("Population after Iteration {}".format(iteration))
    ax.set_xlabel(r"Total Fiber Volume $[\mathrm{\mu m^3}]$")
    ax.set_ylabel("Stress L1-Norm [pPa]")
    loc = "upper right" if limits is not None else "best"
    ax.legend(
        loc=loc,
        handles=[
            mlines.Line2D(
                [0],
                [0],
                marker="o",
                color=mcm.get_cmap("Blues")(0.7),
                label="Selected",
                lw=0,
            ),
            mlines.Line2D(
                [0],
                [0],
                marker="o",
                markeredgecolor=red,
                fillstyle="none",
                label="Removed",
                lw=0,
            ),
            mlines.Line2D([0], [0], marker="x", color="grey", label="Created", lw=0),
        ],
    )

    # Axis limits
    if ax.get_xlim()[0] < 0.0:
        ax.set_xlim(left=0.0)
    if ax.get_ylim()[0] < 0.0:
        ax.set_ylim(bottom=0.0)
    if limits is not None:
        # ax.set_xlim(left=limits["volume"]["min"], right=limits["volume"]["max"])
        # ax.set_ylim(bottom=limits["stress"]["min"], top=limits["stress"]["max"])
        cb.ax.set_ylim(bottom=limits["fitness"]["min"], top=limits["fitness"]["max"])

    fig.savefig(outfile)
    plt.close(fig)


def entrypoint():
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(
        description="DUNE STRUCTURES Population Fitness Plot CLI"
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
        "--outfilename", "-o", help="Output file base name", default="population"
    )
    args = parser.parse_args()

    # Load data
    print("Loading CSV data...")
    results = load_iteration_results(args.root_dir)
    limits = determine_plot_limits(results)

    # Plot population
    print("Plotting populations for {} iterations...".format(len(results)))
    for i, df in enumerate(results):
        plot_population(
            df,
            i,
            os.path.join(
                args.root_dir,
                "{}-{:03d}.{}".format(args.outfilename, i, args.extension),
            ),
            limits,
        )
