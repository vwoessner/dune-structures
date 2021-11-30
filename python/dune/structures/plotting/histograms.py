import os
from matplotlib import markers

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import numpy as np
from cycler import cycler

from ._utils import load_fibers, load_iteration_results
from ._tol_colors import tol_cset


def extract_fiber_data(df_results):
    ret = pd.DataFrame(columns=["length", "radius", "volume"])
    for config in df_results.itertuples():
        fibers = load_fibers(config.filepath)
        length = [
            np.sqrt(np.sum((np.array(fiber["end"]) - np.array(fiber["start"])) ** 2))
            for fiber in fibers
        ]
        radius = [fiber["radius"] for fiber in fibers]
        volume = [np.pi * le * ra ** 2 for le, ra in zip(length, radius)]
        ret = pd.concat(
            [
                ret,
                pd.DataFrame.from_dict(
                    {"length": length, "radius": radius, "volume": volume}
                ),
            ],
            ignore_index=True,
        )

    return ret


def plot_histogram(root_dir, iterations, num_samples, extension="pdf"):
    # Load result data for all iterations
    print("Loading CSV data...")
    results = load_iteration_results(root_dir)

    # Select configuration samples
    print("Loading fiber data from selected YAML config files...")
    data_plot = pd.DataFrame(
        columns=["iteration", "sorting", "length", "radius", "volume"]
    )
    for iter in iterations:
        # Extract fiber data
        results_stress = (
            results[iter].sort_values("stress", ascending=True).head(num_samples)
        )
        fiber_data_stress = extract_fiber_data(results_stress)
        fiber_data_stress["sorting"] = "Lowest stress"

        results_fitness = (
            results[iter].sort_values("fitness", ascending=False).head(num_samples)
        )
        fiber_data_fitness = extract_fiber_data(results_fitness)
        fiber_data_fitness["sorting"] = "Highest fitness"

        # Concatenate datasets
        fiber_data = pd.concat(
            [fiber_data_stress, fiber_data_fitness], ignore_index=True
        )
        fiber_data["iteration"] = iter
        data_plot = pd.concat([data_plot, fiber_data], ignore_index=True)

    # Plot the damn thing!
    print("Creating plots...")
    for y, label in [
        ("length", "Fiber Length"),
        ("radius", "Fiber Radius"),
        ("volume", "Fiber Volume"),
    ]:
        fig, ax = plt.subplots(1, 1, constrained_layout=True, dpi=300)
        sns.violinplot(
            x="iteration",
            y=y,
            hue="sorting",
            data=data_plot,
            palette=tol_cset("bright"),
            # split=True,
            scale="count",
            scale_hue=False,
            inner="box",
            ax=ax,
        )
        ax.set_title("Distributions of " + label)
        ax.set_xlabel("Iteration")
        exp = r"^3" if y == "volume" else r""
        ax.set_ylabel(label + r" $[\mathrm{\mu m}" + exp + r"]$")
        ax.set_ylim(bottom=0.0)
        ax.legend(title=None)
        fig.savefig(os.path.join(root_dir, "dist-{}.{}".format(y, extension)))
        plt.close(fig)

    # Plot a jointplot
    # g = sns.jointplot(
    #     data=data_plot.loc[data_plot["sorting"] == "Highest fitness"],
    #     x="length",
    #     y="radius",
    #     hue="iteration",
    #     palette=list(tol_cset("bright"))[:len(iterations)],
    # )
    # g.plot(sns.scatterplot, sns.histplot)
    # g.figure.savefig(os.path.join(root_dir, "dist-corr.{}".format(extension)))

    fig, ax = plt.subplots(1, 1, constrained_layout=True, dpi=300)
    markersize = mpl.rcParams["lines.markersize"]
    hc = tol_cset("high-contrast")
    styles = cycler(
        color=[hc.blue, hc.red, hc.yellow],
        markersize=np.array([1.3, 1.0, 0.7]) * markersize,
    )
    data_plot = data_plot.loc[data_plot["sorting"] == "Highest fitness"]
    # print(data_plot)
    for iter, style in zip(iterations, styles):
        df = data_plot.loc[data_plot["iteration"] == iter]
        # print(df)
        ax.plot(df["length"], df["radius"] * 1e3, "o", label=str(iter), **style)
    ax.set_xlabel(r"Fiber Lengh $[\mathrm{\mu m}]$")
    ax.set_ylabel(r"Fiber Radius $[\mathrm{n m}]$")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.legend(title="Iteration")
    ax.set_title("Fiber properties for highest fitness")
    fig.savefig(os.path.join(root_dir, "dist-corr.{}".format(extension)))


def entrypoint():
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description="DUNE STRUCTURES Histogram Plot CLI")
    parser.add_argument(
        "root_dir",
        type=os.path.realpath,
        help="Root directory of the optimization algorithm output",
    )
    parser.add_argument("--iterations", "-i", type=int, nargs="*")
    parser.add_argument("--num-samples", "-n", type=int, default=20)
    parser.add_argument(
        "--extension", "-e", help="Filename extension without dot", default="pdf"
    )
    args = parser.parse_args()

    # Call function
    mpl.rcParams["figure.dpi"] = 300
    plot_histogram(args.root_dir, args.iterations, args.num_samples, args.extension)
