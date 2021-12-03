import os
import glob

from ruamel.yaml import YAML
import pandas as pd


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


def load_results(filepath):
    """Load the results of an optimization iteration into a Pandas dataframe"""
    return pd.read_csv(filepath, header=0)


def load_iteration_results(root_dir, results_file="results.csv"):
    return [
        load_results(os.path.join(iter_dir, results_file))
        for iter_dir in sorted(glob.iglob(os.path.join(root_dir, "iteration-*/")))
    ]


def load_fibers(yaml_file_path):
    """Load the fiber information from a run config file"""
    cfg = load_yaml_config(yaml_file_path)
    fiber_block = get_block_by_name_recursive(
        cfg["solver"]["blocks"], "linearsolver_0"
    )["operator"]["reinforced_operator"]["fibres"]
    if fiber_block is not None:
        return [fiber for fiber in fiber_block]
    else:
        return []
