import os
from collections import namedtuple

import numpy as np
from numpy.polynomial import Polynomial
from ruamel.yaml import YAML

import matplotlib.pyplot as plt

from scipy.optimize import minimize
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max

from dune.structures import VTKVertexReader, VTKCellReader
from dune.structures.fibergrowth.line import Line, three_line_integral, angle_between


Fiber = namedtuple(
    "Fiber",
    ["start", "end", "radius", "youngs_modulus", "prestress", "default"],
    defaults=[False],
)


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


def choose_from_local_maxima(
    positions, values, precision, min_distance, filter_sigma, filename
):
    # Create coordinate mesh
    x_min, x_max = np.amin(positions[..., 0]), np.amax(positions[..., 0])
    y_min, y_max = np.amin(positions[..., 1]), np.amax(positions[..., 1])
    x_prec = int(abs(x_max - x_min) * precision)
    y_prec = int(abs(y_max - y_min) * precision)
    x = np.linspace(x_min, x_max, x_prec)
    y = np.linspace(y_min, y_max, y_prec)
    xx, yy = np.meshgrid(x, y)

    # Interpolate data into mesh
    values_grid = griddata(positions, values, (xx, yy), method="linear", fill_value=0.0)

    # Add a gaussian filter
    values_grid = gaussian_filter(values_grid, sigma=filter_sigma * precision)

    # print(values_grid.shape)
    # print(x_prec, precision * min_distance)

    # Evaluate local maxima
    idx = peak_local_max(
        values_grid,
        min_distance=int(precision * min_distance),
        exclude_border=False,
        threshold_rel=0.1,
        p_norm=2,  # Euclidian distance
    )

    plt.imshow(values_grid)
    plt.plot(idx[..., 1], idx[..., 0], "ro")

    name, _ = os.path.splitext(filename)
    filename = name + "-max.pdf"
    plt.savefig(filename)
    plt.close()

    # Return positions and values
    # print(values_grid.shape)
    # print(xx.shape, yy.shape)
    # print(idx.shape)
    # print(idx)
    x_values = xx[idx[..., 0], idx[..., 1]]
    y_values = yy[idx[..., 0], idx[..., 1]]
    values = values_grid[idx[..., 0], idx[..., 1]]
    points = np.array([[x, y] for x, y in zip(x_values, y_values)])
    return points, values


def load_stress_ev(filename, dataset="stress_ev_1_"):
    vtkfile = VTKCellReader(filename)

    return vtkfile.cell_centers, vtkfile[dataset]


def load_stress_mises(filename, dataset="vonmises_"):
    try:
        vtkfile = VTKVertexReader(filename)
        return vtkfile.points, vtkfile[dataset]
    except:
        vtkfile = VTKCellReader(filename)
        return vtkfile.cell_centers, vtkfile[dataset]


def recombine_fibers(fibers, max_distance, max_angle, **kwargs):
    if fibers is None:
        return []
    elif len(fibers) < 2:
        return fibers

    # Degree to rad!
    import math

    max_angle = max_angle * math.pi / 180

    def line_from_entry(idx):
        return Line(fibers[idx]["start"], fibers[idx]["end"])

    # Find nearest neighbors
    neighbors = []
    for i in range(len(fibers)):
        root = line_from_entry(i)
        dist = np.inf
        idx = i + 1

        for j in range(i + 1, len(fibers)):
            fiber = line_from_entry(j)

            def distance(point1, point2):
                return np.sqrt(np.sum((point1 - point2) ** 2))

            def min_distance(line1, line2):
                return min(
                    distance(line1.start, line2.start),
                    distance(line1.end, line2.end),
                    distance(line1.start, line2.end),
                    distance(line1.end, line2.start),
                )

            dist_new = min_distance(root, fiber)
            if dist_new < dist:
                dist = dist_new
                target = fiber
                idx = j

        neighbors.append((i, idx, dist, angle_between(root, target)))

    print("All NNs:")
    print(neighbors)

    # Apply filtering
    neighbors = sorted(
        list(
            filter(
                lambda x: np.isfinite(x[2])
                and x[2] <= max_distance
                and np.isfinite(x[3])
                and x[3] <= max_angle,
                neighbors,
            )
        ),
        key=lambda x: x[2],
    )

    # Return now if no fiber pairs made it to this point
    if not neighbors:
        return fibers

    # Throw out duplicates with longer distances
    index_set = {neighbors[0][0], neighbors[0][1]}
    for pair_info in list(neighbors[1:]):
        if pair_info[0] in index_set or pair_info[1] in index_set:
            neighbors.remove(pair_info)
        index_set.update({pair_info[0], pair_info[1]})

    # Return now if no fiber pairs made it to this point
    if not neighbors:
        return fibers

    print("Processed NNs:")
    print(neighbors)

    # Compute new fibers
    # def func(line_coords, line1, line2):
    #     line_new = Line(
    #         [line_coords[0], line_coords[1]], [line_coords[2], line_coords[3]]
    #     )
    #     return three_line_integral(line_new, line1, line2)

    fibers_to_delete = []
    for pair_info in neighbors:
        idx_1 = pair_info[0]
        idx_2 = pair_info[1]
        line1 = line_from_entry(idx_1)
        line2 = line_from_entry(idx_2)

        # NEED: Polynomial fit here!
        # see https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.Polynomial.fit.html
        x = [line1.start[0], line1.end[0], line2.start[0], line2.end[0]]
        y = [line1.start[1], line1.end[1], line2.start[1], line2.end[1]]
        fit = Polynomial.fit(x=x, y=y, deg=1)

        x_start = np.amin(x)
        x_end = np.amax(x)
        new_line = Line([x_start, fit(x_start)], [x_end, fit(x_end)])

        # Increase radius if fiber length is not increased significantly
        fiber_data = dict(**fibers[idx_1])
        fiber_data_2 = dict(**fibers[idx_2])
        length_init = max(line1.length, line2.length)
        radius = max(fiber_data["radius"], fiber_data_2["radius"])
        if new_line.length < length_init * 1.1:
            radius = fiber_data["radius"] + fiber_data_2["radius"]
        fiber_data["radius"] = radius

        # start_mean = np.mean([line1.start, line2.start], axis=0)
        # end_mean = np.mean([line1.end, line2.end], axis=0)

        # dist = distance(start_mean, end_mean)
        # if dist < line1.length or dist < line2.length:
        #     continue

        # res = minimize(func, x0=[start_mean[0], start_mean[1], end_mean[0], end_mean[1]], args=(line1, line2))
        # pos = res.x
        # new_line = Line([pos[0], pos[1]], [pos[2], pos[3]])
        # new_line = Line([start_mean[0], start_mean[1]], [end_mean[0], end_mean[1]])

        fiber_data["start"] = [new_line.start.item(i) for i in range(2)]
        fiber_data["end"] = [new_line.end.item(i) for i in range(2)]
        new_fiber = Fiber(**fiber_data)
        fibers.append(new_fiber._asdict())
        fibers_to_delete.extend([idx_1, idx_2])

    # Delete old fibers
    # print("Deleting fibers:")
    # print(fibers_to_delete)
    fibers_to_delete = sorted(fibers_to_delete, reverse=True)
    print("Deleting fibers:")
    print(fibers_to_delete)
    for idx in fibers_to_delete:
        fibers.pop(idx)

    return fibers


def compute_new_fibers(
    datafile,
    ev_dataset,
    scale,
    radius,
    youngs_modulus,
    min_length,
    max_length,
    prestress,
    fib_create_dist,
    gauss_sigma,
    **kwargs
):
    # Load the stress EVs
    points, stress_vectors = load_stress_ev(datafile, ev_dataset)

    # Make sure stress EVs are normalized
    stress_vectors = (
        stress_vectors.T / np.sqrt(np.sum(stress_vectors ** 2, axis=-1))
    ).T

    # Compute rectangular bounds
    bounds = np.zeros((3, 2))
    for idx in range(3):
        bounds[idx] = [np.min(points[:, idx]), np.max(points[:, idx])]

    def is_inside(pos):
        for idx in range(3):
            if pos[idx] < bounds[idx, 0] or pos[idx] > bounds[idx, 1]:
                return False
        return True

    # Collapse points into 2D
    points = points[..., 0:2]

    # Choose start positions of new fibers
    # idx_new = choose_from_magnitude(stress_vectors)
    # print("Possible fiber start positions:")
    # print(points[idx_new])

    pp, stress_magnitudes = load_stress_mises(datafile)
    # stress_magnitudes = np.sqrt(np.sum(stress_vectors ** 2, axis=-1))
    pos, values = choose_from_local_maxima(
        pp[..., 0:2],
        stress_magnitudes,
        precision=20,
        min_distance=fib_create_dist,
        filename=datafile,
        filter_sigma=gauss_sigma,
    )

    # Find closest location of VTK values
    idx_new = []
    for point in pos:
        distances = np.sqrt(np.sum((points - point) ** 2, axis=1))
        idx_new.append(np.argmin(distances))

    # Create fibers
    scale = scale
    fibers = []
    for idx, magnitude in zip(idx_new, values):
        center = points[idx]
        direction = scale * magnitude * stress_vectors[idx, 0:2]  # 2D!

        length = np.sqrt(np.sum(direction ** 2))
        if length > max_length:
            direction = direction * max_length / length
        elif length < min_length:
            direction = direction * min_length / length

        start = center - direction / 2
        end = start + direction / 2

        # Outside? Then flip the fiber
        # if not is_inside(end):
        #     end = start - scale

        # Can this happen? Skip it anyway...
        # if not is_inside(end):
        #     continue

        # Make YAML understand the coordinates
        start = [start.item(i) for i in range(2)]
        end = [end.item(i) for i in range(2)]
        # start, end = list(start)[:2], list(end)[:2]
        fibers.append(Fiber(start, end, radius, youngs_modulus, prestress)._asdict())

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


def update_yaml_data(yaml_data, fibers, overwrite=False):
    yaml_op = get_block_by_name(yaml_data, "linearsolver_0")["operator"][
        "reinforced_operator"
    ]
    if overwrite or "fibres" not in yaml_op or yaml_op["fibres"] is None:
        yaml_op["fibres"] = fibers
    else:
        yaml_op["fibres"].extend(fibers)
    return yaml_data


def get_yaml_data_fibers(yaml_data):
    yaml_op = get_block_by_name(yaml_data, "linearsolver_0")["operator"][
        "reinforced_operator"
    ]
    return yaml_op["fibres"]


def load_fibergrowth_data(yaml_data):
    return yaml_data["fibergrowth"]


def format_iteration(string, iteration):
    return string.format(iteration="{:03d}".format(iteration))


def fibergrowth_cfg(yaml_input_file, create=True, recombine=True):
    """Run the fiber growth algorithm from a YAML input file which is assumed to have
    been executed by a simulation, i.e., with existing output.
    """

    # Load the input file and data
    yaml_data = load_yaml_config(yaml_input_file)
    fibergrowth_data = load_fibergrowth_data(yaml_data)
    vtk_output = get_block_by_name(yaml_data, "visualization_0")["filename"] + ".vtu"

    if recombine:
        # Fetch all current fibers from YAML file
        fiber_data = get_yaml_data_fibers(yaml_data)
        fibers = recombine_fibers(fiber_data, **fibergrowth_data["fibers"])
        yaml_data = update_yaml_data(yaml_data, fibers, overwrite=True)

    # Compute new fibers
    if create:
        fibers = compute_new_fibers(vtk_output, **fibergrowth_data["fibers"])
        yaml_data = update_yaml_data(yaml_data, fibers)

    return yaml_data
