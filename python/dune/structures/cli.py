import os
import argparse
import subprocess
import logging

from dune.structures.fibergrowth.evaluate import (
    load_yaml_config,
    load_fibergrowth_data,
    format_iteration,
    write_yaml_config,
    fibergrowth_cfg,
    get_block_by_name,
)

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
