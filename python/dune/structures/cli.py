import os
import argparse
import subprocess
import logging

# TODO: Hardcoded path! Should be configured by CMake instead
APP_PATH = "{}/../../../build-cmake/apps".format(os.path.dirname(__file__))


def run(args, logger):
    # Check the executable
    exe = os.path.abspath(os.path.join(APP_PATH, args.executable))
    if not os.path.isfile(exe):
        raise FileNotFoundError("Executable not found: {}".format(exe))

    # Check the input file
    input_file = os.path.abspath(args.input_file)
    if not os.path.isfile(input_file):
        raise FileNotFoundError("Input file not found: {}".format(input_file))

    # Execute the process, capturing the error pipe
    process_args = [exe, "run", input_file]
    logger.debug("Executing process: {}".format(process_args))
    output = subprocess.check_output(process_args, stderr=subprocess.STDOUT, text=True)

    # Report the process output
    log_file = args.log_file
    if log_file is not None:
        with open(log_file, "w") as file:
            f.write(output)
    else:
        logger.info("Process output:")
        print(output)


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

    return parser


def cli():
    # Retrieve the logger
    logger = logging.getLogger("CLI")
    logger.setLevel(logging.INFO)

    # Parse arguments
    parser = create_argparse()
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level.upper())
    if args.action == "run":
        run(args, logger)
    else:
        logger.warning("No command line arguments given!")
        parser.print_help()
