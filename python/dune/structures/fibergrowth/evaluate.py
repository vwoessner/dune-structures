from collections import namedtuple

import numpy as np
from ruamel.yaml import YAML

from dune.structures import VTKVertexReader


def get_block_by_name(yaml_data, blockname):
    for block in yaml_data["solver"]["blocks"]:
        if block["_blockname"] == blockname:
            return block


def get_blocks_by_type(yaml_data, blocktype):
    Block = namedtuple(Block, ["name", "data"])
    return [
        Block(block["_blockname"], block)
        for block in yaml_data["solver"]["blocks"]
        if block["_type"] == blocktype
    ]


def choose_from_magnitude(
    eigenvectors, percentage=0.01, threshold=None, *, threshold_type="absolute"
):
    """Takes a list of eigenvectors and returns the indices of the ones with the
    highest magnitude, according to the input queries.
    """
    # Order by descending magnitude
    magnitude = np.sqrt(np.sum(eigenvectors ** 2, axis=-1))
    idx = np.flip(np.argsort(magnitude))

    # Choose top percentage
    return idx[: int(len(idx) * percentage)]


def load_data(filename, dataset="stress_ev_1_"):
    vtkfile = VTKVertexReader(filename)

    return vtkfile.points, vtkfile[dataset]


def compute_new_fibers(datafile, ev_dataset, scale, radius, youngs_modulus, min_length, max_length):
    # Load the data
    points, stress_vectors = load_data(datafile, ev_dataset)

    # Compute rectangular bounds
    bounds = np.zeros((3, 2))
    for idx in range(3):
        bounds[idx] = [np.min(points[:, idx]), np.max(points[:, idx])]

    def is_inside(pos):
        for idx in range(3):
            if pos[idx] < bounds[idx, 0] or pos[idx] > bounds[idx, 1]:
                return False
        return True

    # Choose start positions of new fibers
    idx_new = choose_from_magnitude(stress_vectors)
    print("Possible fiber start positions:")
    print(points[idx_new])

    # Create fibers
    Fiber = namedtuple("Fiber", ["start", "end", "radius", "youngs_modulus"])
    scale = scale
    fibers = []
    for idx in idx_new:
        start = points[idx]
        direction = scale * stress_vectors[idx]
        end = start + direction

        length = np.sqrt(np.sum((end - start) ** 2))
        if length > max_length:
            direction = direction * max_length / length
        elif length < min_length:
            direction = direction * min_length / length

        # Outside? Then flip the fiber
        if not is_inside(end):
            end = start - scale

        # Can this happen? Skip it anyway...
        if not is_inside(end):
            continue

        # Make YAML understand the coordinates
        start = [start.item(i) for i in range(2)]
        end = [end.item(i) for i in range(2)]
        # start, end = list(start)[:2], list(end)[:2]
        fibers.append(Fiber(start, end, radius, youngs_modulus)._asdict())

    # print("New fibers:")
    # for fiber in fibers:
        # print("Start: {}, End: {}".format(fiber["start"], fiber["end"]))
    return fibers


def load_yaml_config(filename):
    with open(filename, "r") as file:
        yaml = YAML(typ="safe")
        yaml.allow_duplicate_keys = True
        return yaml.load(file)


def write_yaml_config(filename, yaml_data):
    with open(filename, "w") as file:
        yaml = YAML(typ="safe")
        yaml.allow_duplicate_keys = True
        yaml.default_flow_style = False
        yaml.dump(yaml_data, file)


def update_yaml_data(yaml_data, fibers):
    yaml_op = get_block_by_name(yaml_data, "linearsolver_0")["operator"][
        "reinforced_operator"
    ]
    if not "fibres" in yaml_op or yaml_op["fibres"] is None:
        yaml_op["fibres"] = fibers
    else:
        yaml_op["fibres"].extend(fibers)
    return yaml_data


def load_fibergrowth_data(yaml_data):
    return yaml_data["fibergrowth"]


def format_iteration(string, iteration):
    return string.format(iteration="{:03d}".format(iteration))


def fibergrowth_cfg(yaml_input_file):
    """Run the fiber growth algorithm from a YAML input file which is assumed to have
    been executed by a simulation, i.e., with existing output.
    """

    # Load the input file and data
    yaml_data = load_yaml_config(yaml_input_file)
    fibergrowth_data = load_fibergrowth_data(yaml_data)
    vtk_output = get_block_by_name(yaml_data, "visualization_0")["filename"] + ".vtu"

    # Evaluate solution and compute new fibers
    fibers = compute_new_fibers(vtk_output, **fibergrowth_data["fibers"])
    return update_yaml_data(yaml_data, fibers)
