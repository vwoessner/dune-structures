import re
import numpy as np
from numpy.random import default_rng, shuffle
import subprocess as sp
import os
import copy
import multiprocessing as mp

from ruamel.yaml import YAML

from .line import Line
from .evaluate import (
    Fiber,
    get_block_by_name,
    write_yaml_config,
    load_yaml_config,
    load_stress_ev,
    load_stress_mises,
    choose_from_local_maxima,
)

# OptData = namedtuple(
#     "OptData", ["population_size", "x_bounds", "y_bounds", "max_fibers"]
# )

REGEX = r"stress l2: (-?[\d.]+e(?:\+|\-)?\d+)"


def transpose_bounds(x_bounds, y_bounds):
    bounds = np.vstack((x_bounds, y_bounds)).T
    return bounds[0], bounds[1]


def random_fiber(x_bounds, y_bounds):
    low, high = transpose_bounds(x_bounds, y_bounds)
    rng = default_rng()
    start, end = rng.uniform(low, high, size=(2, 2))
    return Line(start, end)


def random_fiber_from_stress(
    yaml_config_file, gauss_sigma, scale, fib_create_dist, rng
):
    # Load VTK
    directory = os.path.dirname(yaml_config_file)
    cfg = load_yaml_config(yaml_config_file)
    vtk_output = os.path.join(
        directory, get_block_by_name(cfg, "visualization_0")["filename"] + ".vtu"
    )

    # Load stresses
    points, stress_ev = load_stress_ev(vtk_output)
    points_mises, stress_mises = load_stress_mises(vtk_output)

    # Make sure stress EVs are normalized
    stress_ev = (stress_ev.T / np.sqrt(np.sum(stress_ev ** 2, axis=-1))).T

    # Collapse points into 2D
    points = points[..., 0:2]

    # Evaluate local maxima
    pos, values = choose_from_local_maxima(
        points_mises[..., 0:2],
        stress_mises,
        precision=20,
        min_distance=fib_create_dist,
        filename=vtk_output,
        filter_sigma=gauss_sigma,
    )

    # Find closest location of VTK values
    idx_new = []
    for point in pos:
        distances = np.sqrt(np.sum((points - point) ** 2, axis=1))
        idx_new.append(np.argmin(distances))

    # Random selection with magnitude as weight
    idx = rng.choice(
        list(range(len(idx_new))),
        p=np.log(values) / np.sum(np.log(values)),  # LOG for scaling!
        replace=False,
        shuffle=False,
    )

    # Set Fiber
    center = points[idx_new[idx]]
    length = scale * values[idx]
    direction = length * stress_ev[idx_new[idx], 0:2]  # 2D!
    return Line(center + direction / 2, center - direction / 2)


def load_optimization_data(filename):
    yaml = YAML(typ="safe")
    with open(filename, "r") as file:
        return yaml.load(file)


def generate_initial_population(data, rng):
    num_fibers = rng.integers(
        data["min_num_fibers"], data["max_num_fibers"], size=data["population_size"]
    )
    population = [
        [random_fiber(data["x_bounds"], data["y_bounds"]) for i in range(num)]
        for num in num_fibers
    ]
    return population


def crossover_single(parent1, parent2, rng):
    rng = default_rng()
    cross_point = rng.integers(1, min(len(parent1), len(parent2)))
    # print("Crossover point: {}".format(cross_point))
    child1 = list(parent1[:cross_point]) + parent2[cross_point:]
    child2 = list(parent2[:cross_point]) + parent1[cross_point:]
    # rng.shuffle(child1)
    # rng.shuffle(child2)
    return [child1, child2]


def crossover(population, data, rng):
    new_pop = []
    while len(new_pop) < data["population_size"]:
        select = rng.integers(0, len(population), size=2)
        new_pop.extend(
            crossover_single(population[select[0]], population[select[1]], rng)
        )
    return new_pop


def mutation(population, yaml_file_paths, data, rng):
    for genome, yaml_file in zip(population, yaml_file_paths):
        rand = rng.uniform(size=2)
        if (
            rand[0] < data["fiber_create_probability"]
            and len(genome) < data["max_num_fibers"]
        ):
            # Create fiber
            # genome.append(random_fiber(data["x_bounds"], data["y_bounds"]))
            genome.append(
                random_fiber_from_stress(
                    yaml_file,
                    data["gauss_sigma"],
                    data["scale"],
                    data["fib_create_dist"],
                    rng,
                )
            )
        if (
            rand[1] < data["fiber_delete_probability"]
            and len(genome) > data["min_num_fibers"]
        ):
            # Delete random fiber
            genome.pop(rng.integers(0, len(genome)))

        # Change fibers
        rand = rng.uniform(size=len(genome))
        # print(rand)
        var = data["fiber_mutation_stddev"] ** 2
        cov = [[var, 0], [0, var]]
        prob = data["fiber_mutate_probability"]
        # print(rand < prob)
        # print((rand < prob).nonzero())
        to_mutate = rand < prob
        # print(to_mutate)
        if to_mutate.any():
            for idx in np.argwhere(to_mutate)[0]:
                gene = genome[idx]
                start = gene.start
                start = rng.multivariate_normal(start, cov)
                end = gene.end
                end = rng.multivariate_normal(end, cov)
                genome[idx] = Line(start, end)


def selection(population, scores, data, rng):
    def pareto_front_indices(scores):
        """Return a boolean array selecting scores in the pareto front"""
        idx = np.full(len(scores), True)
        for i in range(len(scores)):
            for j in range(len(scores)):
                if all(scores[j] <= scores[i]) and any(scores[j] < scores[i]):
                    idx[i] = False
                    break
        return idx

    def crowding_distances(scores):
        # print(scores.shape)
        normed_scores = scores - scores.min(axis=0) / scores.ptp(axis=0)
        distances = np.zeros_like(scores)
        pop_size = scores.shape[0]
        for i in range(scores.shape[1]):
            target_scores = normed_scores[:, i]
            idx = np.argsort(target_scores)
            target_scores = target_scores[idx]
            crowding = np.ones(pop_size)
            crowding[1 : pop_size - 1] = (
                target_scores[2:pop_size] - target_scores[0 : pop_size - 2]
            )
            idx_rev = np.argsort(idx)
            distances[:, i] = crowding[idx_rev]
        return distances

    # Select Pareto fronts until selection size is reached
    # print(scores.shape)
    selection_size = int(data["selection_ratio"] * len(population))
    selected_idx = np.full(len(scores), False)
    while np.count_nonzero(selected_idx) < selection_size:
        # print(selected_idx.shape)
        reverse_idx = np.argwhere(~selected_idx)
        # print(reverse_idx.shape)
        last_pareto = pareto_front_indices(scores[~selected_idx])
        # print(last_pareto.shape)
        selected_idx[reverse_idx[last_pareto]] = True
        # print(selected_idx)

    # print(reverse_idx[last_pareto])
    # print(scores[reverse_idx[last_pareto]])
    # Throw out elements of last front according to crowding distance
    # crowding_mask = np.full_like(selected_idx, True)
    # pop_size = np.count_nonzero(selected_idx)
    deletable_idx = np.squeeze(reverse_idx[last_pareto])
    while np.count_nonzero(selected_idx) > selection_size:
        # Competition between two random individuals
        crowding = crowding_distances(scores[deletable_idx])
        select = rng.integers(0, len(crowding), size=2)
        last_pareto_idx = np.argwhere(last_pareto)
        if np.sum(crowding[select[0]]) > np.sum(crowding[select[1]]):
            # selected_idx[last_pareto_idx[select[1]]] = False
            # last_pareto[select[1]] = False
            delete_idx = select[1]
        else:
            # selected_idx[last_pareto_idx[select[0]]] = False
            # last_pareto[select[0]] = False
            delete_idx = select[0]
        # print(delete_idx)
        selected_idx[deletable_idx[delete_idx]] = False
        deletable_idx = np.delete(deletable_idx, delete_idx)

        # reverse_idx = np.delete(reverse_idx, select)
        # last_pareto = np.delete(last_pareto, select_delete)
        # last_pareto[last_pareto_idx[select_delete]] = False
        # crowding_mask[reverse_idx[last_pareto]] = False
        # selected_idx[reverse_idx[select_delete]] = False

    # selected_idx = selected_idx & crowding_mask
    # print(selected_idx)

    # return population[selected_idx]
    return [population[idx] for idx in selected_idx if idx]

    # idx = np.argsort(stresses)
    # print("Best genomes:")
    # print(idx[: int(data["selection_ratio"] * len(population))])
    # return [
    #     population[i] for i in idx[: int(data["selection_ratio"] * len(population))]
    # ]


def evaluate_fiber_lengths(population, data):
    return [np.sum([fiber.length for fiber in genome]) for genome in population]


def run_and_evaluate_parallel(executable, yaml_file_paths, processes=1):
    # def run_eval(path):
    #     return run_and_evaluate(executable, path)

    with mp.Pool(processes) as pool:
        return pool.starmap(
            run_and_evaluate, [(executable, path) for path in yaml_file_paths]
        )


def run_and_evaluate(executable, yaml_file_path):
    directory = os.path.dirname(yaml_file_path)
    output = sp.check_output(
        [executable, "run", yaml_file_path], cwd=directory, encoding="utf8"
    )
    match = re.search(REGEX, output)
    # print(output, match.group(1))
    return float(match.group(1))


def genome_to_cfg(yaml_cfg, genome):
    fibers = [
        Fiber(
            start=[fiber.start.item(i) for i in range(2)],
            end=[fiber.end.item(i) for i in range(2)],
            **yaml_cfg["genetic_optimization"]["fibers"]
        )._asdict()
        for fiber in genome
    ]
    get_block_by_name(yaml_cfg, "linearsolver_0")["operator"]["reinforced_operator"][
        "fibres"
    ] = fibers


def write_population_cfgs(directory, filename_base, default_yaml_cfg, population):
    if not os.path.isdir(directory):
        os.makedirs(directory)

    filename, ext = os.path.splitext(filename_base)
    file_paths = []
    for idx, genome in enumerate(population):
        yaml_cfg = copy.deepcopy(default_yaml_cfg)
        genome_to_cfg(yaml_cfg, genome)
        outfile = str(get_block_by_name(yaml_cfg, "visualization_0")["filename"])
        get_block_by_name(yaml_cfg, "visualization_0")["filename"] = outfile.format(idx)
        filepath = os.path.join(directory, filename + "-{:03d}".format(idx) + ext)
        write_yaml_config(filepath, yaml_cfg)
        file_paths.append(filepath)

    return file_paths


# def load_optimization_data(yaml_data):
#     return yaml_data["genetic_optimization"]


# def random_fiber(x_bounds, y_bounds):
#     bounds = np.vstack((x_bounds, y_bounds)).T
#     low, high = bounds[0], bounds[1]

#     rng = default_rng()
#     coords = rng.uniform(low=low, high=high, size=(2, 2))

#     return Line(coords[0], coords[1])


# def initialize_population(opt_data):
#     # Determine number of fibers for each pop
#     rng = default_rng()
#     num_fibers = rng.integers(
#         low=1, high=opt_data.max_fibers, endpoint=True, size=opt_data.population_size
#     )
#     population = [
#         [random_fiber(opt_data.x_bounds, opt_data.y_bounds) for i in range(num)]
#         for num in num_fibers
#     ]
