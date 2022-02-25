import os
import glob
import argparse
import subprocess
import logging
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from numpy.random import default_rng
import pandas as pd

from dune.structures.fibergrowth.evaluate import (
    load_yaml_config,
    load_fibergrowth_data,
    format_iteration,
    write_yaml_config,
    fibergrowth_cfg,
    get_block_by_name,
)

from dune.structures.fibergrowth.genetic_opt import (
    fill_random_population,
    write_population_cfgs,
    run_and_evaluate_parallel,
    selection,
    crossover,
    mutation,
    evaluate_fiber_lengths,
    evaluate_fiber_volumes,
    get_mesh_triangulation,
)

from dune.structures.plotting.displacement import plot_mean_displacement
from dune.structures.plotting.population import plot_population

# TODO: Hardcoded path! Should be configured by CMake instead
APP_PATH = "{}/../../../build-cmake/apps".format(os.path.dirname(__file__))


def run(executable, input_file, logger, log_file=None, **kwargs):
    # Check the executable
    exe = os.path.abspath(os.path.join(APP_PATH, executable))
    if not os.path.isfile(exe):
        raise FileNotFoundError("Executable not found: {}".format(exe))

    # Check the input file
    input_file = os.path.abspath(input_file)
    if not os.path.isfile(input_file):
        raise FileNotFoundError("Input file not found: {}".format(input_file))

    # Execute the process, capturing the error pipe
    process_args = [exe, "run", input_file]
    logger.debug("Executing process: {}".format(process_args))
    try:
        output = subprocess.check_output(
            process_args, stderr=subprocess.STDOUT, text=True
        )
    except subprocess.CalledProcessError as exc:
        output = exc.output
        logger.error("Executing the process failed:")
        print(output)

    # Report the process output
    if log_file is not None:
        with open(log_file, "w") as file:
            f.write(output)
    else:
        logger.info("Process output:")
        print(output)


def fibergrowth(executable, input_file, logger, **kwargs):
    # Determine file paths
    input_file = os.path.abspath(input_file)
    input_dir, input_filename = os.path.split(input_file)

    # Load the overall input file
    yaml_input = load_yaml_config(input_file)
    fibergrowth_data = load_fibergrowth_data(yaml_input)

    # Write the initial simulation file
    logger.info("Writing initial simulation file")
    output_placeholder = get_block_by_name(yaml_input, "visualization_0")["filename"]
    get_block_by_name(yaml_input, "visualization_0")["filename"] = format_iteration(
        str(output_placeholder), 0
    )
    input_filename, ext = os.path.splitext(input_filename)
    simulation_input_file = os.path.join(input_dir, input_filename + "-000" + ext)
    write_yaml_config(simulation_input_file, yaml_input)

    # Run initial simulation
    logger.info("Executing initial simulation")
    run(executable, simulation_input_file, logger=logger)

    # Run the iterative algorithm
    for iteration in range(1, fibergrowth_data["iterations"] + 1):

        # Compute and write new fibers
        logger.info("Running fibergrowth algorithm iteration {}".format(iteration))

        # Odd
        if iteration % 2:
            yaml_data = fibergrowth_cfg(simulation_input_file, recombine=False)
        # Even
        else:
            yaml_data = fibergrowth_cfg(simulation_input_file, create=False)

        get_block_by_name(yaml_data, "visualization_0")["filename"] = format_iteration(
            str(output_placeholder), iteration
        )
        simulation_input_file = os.path.join(
            input_dir, input_filename + "-{:03d}".format(iteration) + ext
        )

        write_yaml_config(simulation_input_file, yaml_data)

        # Run simulation
        logger.info("Executing simulation {}".format(iteration))
        run(executable, simulation_input_file, logger=logger)


def genetic_opt(executable, input_file, logger, **kwargs):
    # Check the executable
    exe = os.path.abspath(os.path.join(APP_PATH, executable))
    if not os.path.isfile(exe):
        raise FileNotFoundError("Executable not found: {}".format(exe))

    # Load optimization data

    # Determine file paths
    input_file = os.path.abspath(input_file)
    input_dir, input_filename = os.path.split(input_file)

    # Load the overall input file
    yaml_input = load_yaml_config(input_file)
    optimization_data = yaml_input["genetic_optimization"]
    # output_placeholder = get_block_by_name(yaml_input, "visualization_0")["filename"]

    # Fetch the RNG
    rng = default_rng(optimization_data["seed"])

    # Load the mesh triangulation (for position information)
    tri, media = get_mesh_triangulation(
        yaml_input["solver"]["grid"]["filename"].removeprefix("../")
    )
    trifinder = tri.get_trifinder()

    # Generate initial population
    population = []
    fill_random_population(population, optimization_data, rng, trifinder, media)
    population_old = None
    scores_old = None

    # List of dataframes for each iteration
    data_iteration = []

    for it in range(optimization_data["iterations"]):
        print("Iteration {}".format(it))
        # Write YAML input files
        iter_dir = os.path.join(input_dir, "iteration-{:03d}".format(it))
        yaml_file_paths = write_population_cfgs(
            iter_dir, input_filename, yaml_input, population
        )

        # Run and evaluate
        stresses = run_and_evaluate_parallel(
            exe, yaml_file_paths, optimization_data, processes=8
        )
        lengths = evaluate_fiber_volumes(population, optimization_data)
        scores = np.column_stack((stresses, lengths))

        # Store entire population for this round
        data = pd.DataFrame(
            {
                "filepath": [os.path.abspath(path) for path in yaml_file_paths],
                "iteration": it,
                "stress": stresses,
                "length": lengths,
            }
        )

        # Fetch selected population from previous round for egalitarian selection
        if it > 0:
            data_prev = data_iteration[-1]
            data = pd.concat([data, data_prev.loc[data_prev["selected"]]])

        # Selection
        pop_to_select = population
        scores_to_select = scores
        if population_old is not None:
            pop_to_select = population + population_old
            scores_to_select = np.concatenate((scores, scores_old))
        population, fitness, select = selection(
            pop_to_select, scores_to_select, optimization_data, rng
        )
        data["fitness"] = fitness
        data["selected"] = select
        fitness = fitness[select]

        # Write the data to a file
        with open(os.path.join(iter_dir, "results.csv"), "w", newline="") as file:
            data_sorted = data.sort_values(
                by=["fitness", "stress"], ascending=[False, True]
            )
            data_sorted.to_csv(file)

        # Retain population before crossover/mutation as "old"
        population_old = population
        scores_old = scores_to_select[select]

        # Crossover, Mutation
        population = crossover(population, fitness, optimization_data, rng)
        # print(population)
        # NOTE: Population was already changed, stress is outdated!
        # NOTE: Maybe use mean stress?
        mutation(population, yaml_file_paths, optimization_data, rng, trifinder, media)
        fill_random_population(population, optimization_data, rng, trifinder, media)

        # Shuffle genomes
        for genome in population:
            rng.shuffle(genome)

        # Retain iteration data
        data_iteration.append(data)

        # Plot population
        plot_population(
            data, it, os.path.join(input_dir, "population-{:03d}.pdf".format(it))
        )

        # Plot the mean of the 10% of best files
        for data_plot, outname in [
            (data_sorted, "fittest"),
            (data.sort_values("stress", ascending=True), "stress-low"),
            (data.sort_values("length", ascending=True), "fiber-low"),
        ]:
            data_plot = data_plot.iloc[: int(data_plot.shape[0] // 10)]
            plot_yaml_paths = data_plot["filepath"]
            plot_vtk_paths = [
                os.path.splitext(path)[0] + ".vtu" for path in plot_yaml_paths
            ]
            plot_mean_displacement(
                plot_vtk_paths,
                plot_yaml_paths,
                "{}-{:03d}.pdf".format(outname, it),
            )

        # Clean up VTK files of removed genes (only from this iteration)
        if optimization_data["remove_vtk_files"]:
            for file in data["filepath"].loc[
                ~data["selected"] & (data["iteration"] == it)
            ]:
                os.remove(os.path.splitext(file)[0] + ".vtu")

    # print("Deleting VTK output files...")
    # for file in glob.iglob(os.path.join(input_dir, "*/*.vtu")):
    #     os.remove(file)

    # Iterate:
    # Run simulations (parallel)
    # Evaluate fitness (parallel)
    # Selection according to fitness
    # Crossover breeding until population is full
    # Mutation: Vary fibers
    # Mutation: New fibers
    # Write YAML input files for next iteration


def create_argparse():
    # Parent parser
    parser = argparse.ArgumentParser(
        description="DUNE STRUCTURES Command Line Interface"
    )
    parser.add_argument(
        "--log-level",
        choices=["debug", "info", "warning", "error", "critical"],
        default="info",
        help="Logging level of the command line interface",
    )

    # Subparser action object
    sub = parser.add_subparsers(dest="action")

    # "Run" parser
    parse_run = sub.add_parser("run", help="Execute a simulation")
    parse_run.add_argument("executable", type=str, help="Name of the app to execute")
    parse_run.add_argument(
        "input_file",
        type=os.path.realpath,
        help="YAML configuration input file for the simulation",
    )
    parse_run.add_argument(
        "--log-file",
        "-l",
        type=os.path.realpath,
        help="File to write simulation output into",
    )
    parse_run.set_defaults(func=run)

    # "Fibergrowth" parser
    parse_fib = sub.add_parser(
        "fibergrowth", help="Run the iterative fibergrowth algorithm"
    )
    parse_fib.add_argument("executable", type=str, help="Name of the app to execute")
    parse_fib.add_argument(
        "input_file",
        type=os.path.realpath,
        help="YAML configuration input file for the fibergrowth algorithm (META-FILE)",
    )
    parse_fib.set_defaults(func=fibergrowth)

    # "Genetic Optimization" parser
    parse_fib = sub.add_parser(
        "optimization", help="Run the iterative genetic optimization algorithm"
    )
    parse_fib.add_argument("executable", type=str, help="Name of the app to execute")
    parse_fib.add_argument(
        "input_file",
        type=os.path.realpath,
        help="YAML configuration input file for the optimization algorithm (META-FILE)",
    )
    parse_fib.set_defaults(func=genetic_opt)

    return parser


def cli():
    # Retrieve the logger
    logger = logging.getLogger("CLI")
    logger.setLevel(logging.INFO)

    # Parse arguments
    parser = create_argparse()
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level.upper())

    # if args.action == "run":
    #     run(**vars(args), logger=logger)
    # elif args.action == "fibergrowth":
    #     fibergrowth(**vars(args), logger=logger)
    # else:
    #     logger.warning("No command line arguments given!")
    #     parser.print_help()

    # Run the default command
    args.func(**vars(args), logger=logger)
