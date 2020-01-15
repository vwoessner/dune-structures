This module contains the source code for the simulation of single cell mechanics and stress fibre formation executed in the STRUCTURES cluster of excellence at Heidelberg university. For questions write to [dominic.kempf@iwr.uni-heidelberg.de](mailto:dominic.kempf@iwr.uni-heidelberg.de).

# 	Installing the software

There is two principal ways to use this software: As dockerized applications or by building it from source. The former is recommended if you intend to *run* simulations, the latter is necessary if you want to *develop* your own simulators.

## Dockerized setup

Make sure that you have Docker installed on your machine.
You can do so e.g. by following [these instructions from Docker](https://docs.docker.com/install/linux/docker-ce/debian/#install-using-the-repository). Having done so, you can run the container interactively by doing:

```
docker run -ti registry.dune-project.org/dominic/dune-structures
```

More information that should be provided here:
* Entrypoints of the Docker container
* Mounting volumes into the container to retrieve visualization data

## Full installation

This is not an installation guide, but merely documents the specific requirements of dune-structures.
It assumes that you know how to build a Dune software stack.

dune-structures is a Dune module and as such has some external dependencies:

* CMake >= 3.1
* C++17-compliant compiler e.g. gcc >= 8.1
* MPI (e.g. Debian packages `libopenmpi-dev` and `openmpi-bin`)

dune-structures requires the following Dune modules:

* All the [Dune core modules](https://gitlab.dune-project.org/core)
* [dune-typetree](https://gitlab.dune-project.org/staging/dune-typetree.git)
* [dune-functions](https://gitlab.dune-project.org/staging/dune-functions.git)
* [dune-pdelab](https://gitlab.dune-project.org/pdelab/dune-pdelab.git)
* [dune-testtools](https://gitlab.dune-project.org/quality/dune-testtools.git)
* [dune-uggrid](https://gitlab.dune-project.org/staging/dune-uggrid.git)
* [dune-codegen](https://gitlab.dune-project.org/extensions/dune-codegen.git)

Note that dune-codegen requires special care in setup. It needs to
* be cloned with the `--recursive` flag
* You need to run `./patches/apply_patches.sh` from the module root directory after cloning

Currently, the latest versions (master) of all Dune modules are necessary.
Furthermore, dune-structures requires the following external packages:

* A Python interpreter >= 3.6 (e.g. Debian packages `python3-dev`)
* pip (e.g. the Debian package `python3-pip`)
* MuParser (e.g. the Debian package `libmuparser-dev`)
* Gmsh (e.g. Debian packages `gmsh`)

You should build your Dune stack with an options file that contains at least the following entries:

```
CMAKE_FLAGS+="-DDUNE_PYTHON_VIRTUALENV_SETUP=1
              -DDUNE_PYTHON_ALLOW_GET_PIP=1"
```

# Mathematical Abstractions

dune-structures provides solvers for elastostatic and elastodynamic simulations with a variety of linear and nonlinear material laws.

## Basic equations

Given a reference configuration $`T\subseteq\mathbb{R}^3`$ of an elastic body, we aim to calculate the displacement field $`\mathbf{u}:T\rightarrow\mathbb{R}^3`$ that describes the deformation of the body subject to the given forces.

## Strain Measures

Material laws in hyperelasticity are often expressed using a variety of different strain measures.
These are summarized here with the goal to document the implemented physics:

* *Deformation gradient*: $`\mathbf{F} = \mathbb{I} + \nabla\mathbf{u}`$
* *Infinitesimal strain*: $`\frac{1}{2}(\nabla\mathbf{u} + \nabla\mathbf{u}^T)`$
* *Right Cauchy-Green tensor*: $`\mathbf{C} = \mathbf{F}^T\mathbf{F}`$
* *Left Cauchy-Green tensor*: $`\mathbf{B} = \mathbf{F}\mathbf{F}^T`$
* *Cauchy-Green strain tensor*: $`\frac{1}{2}(\mathbf{C}-\mathbb{I})`$
* The *Cauchy-Green invariants*:
    * $`I_1 = tr(\mathbf{B})`$
    * $`I_2 = \frac{1}{2}(tr(\mathbf{B})^2 - tr(\mathbf{B}\mathbf{B})`$
    * $`I_3 = \det \mathbf{B}`$
* The *isochoric Cauchy-Green invariants*:
    * $`\bar{I}_1 = J^{-\frac{2}{3}}I_1`$
    * $`\bar{I}_2 = J^{-\frac{4}{3}}I_2`$
    * $`J=\sqrt{I_3}`$

## Finite Element formulation

In its simplest form, the weak formulation is given as

$`r(v) = \int_\Omega P(u) : v dx`$

withe $`P(u)`$ being the first Piola-Kirchhoff stress tensor.

# Software Design

The software for dune-structures is designed such that it can handle a large variety of problems and numerical techniques without recompilation.
This is done in order to be able to focus on productivity with experiments.

# Solver Component Documentation

This section gives a comprehensive overview of implemented solver components and the configuration values they accept.

## Fundamental simulation building blocks

These solver components handle fundamental building blocks of PDELab simulations.

### Interpolation

This step is registered as `interpolation`. It implements interpolation of a given function into the solution vector. This is necessary for initial conditions of instationary problems as well as in order to implement Dirichlet boundary conditions. In the latter case an extension of the Dirichlet boundary condition needs to be interpolated before the solver is applied.

The interpolation step accepts the following keys:
* `functions` expects a comma separated list of expressions. These expressions are interpolated into the leaves of the function space tree in the given order. For elasticity problems this means that three functions for the displacements in $`x`$, $`y`$ and $`z`$ direction need to be given. Alternatively, a single function can be given, which will then be applied to all leaves of the function tree. Within these functions the strings $`x`$, $`y`$ and $`z`$ can be used to refer to the evaluation point in global coordinates.

The following examples implement interpolation of initial conditions $`\mathbf{u}(\mathbf{x})=0`$ and $`\mathbf{u}(\mathbf{x})=\mathbf{x}`$:

```
[interpolation]
functions = 0.0
```

```
[interpolation]
functions = x, y, z
```

### Constraints Assembly

This step is registered as `constraints`. It implements the assembly of the constraints container. This is e.g. necessary for implementation of Dirichlet boundary conditions.

The constraints assembly step accepts the following keys:
* `functions` expects a comma separated list of expressions as previously described for interpolation. The return value of these functions is a bolean value, where $`0`$ and $`1`$ can be used as literals. A true value results in a Dirichlet-constrained degree of freedom.

The following example implements Dirichlet conditions for the $`z`$ component of the displacement in the $`x-y`$-plane and Neumann boundary conditions elsewhere:

```
[constraints]
functions = 0, 0, z < 1e-8
```

### Time-stepping Loop

This step is registered as `timeloop`.
It controls the basic time stepping loop and introduces parameters `"time"` and `"timestep"` into the parameter system.
This step has substeps, which can be added to through the `add` method.

This timestepping step accepts the following keys:
* `steps` is a list of substeps that will be constructed and added to the time stepper.
* `Tstart` is the starting time
* `Tend` is the stopping time
* `timestep` is the time step

## Elasticity-specific building blocks

This section describes some solver steps that are very specific to the applications of dune-structures.
Most other solver parts could be reused in general PDE applications.

### Nonlinear elasticity solver

This step is registered as `elasticity`.
It applies a Newton solver for the hyperelastic problem defined by the generated operator (see above).
It currently does not accept configuration parameters, although there would be plenty of solver related knobs to expose.

### Material parser

This step is registered as `material`.
It introduces the material class as `"material"` into the parameter system, which is of type `std::shared_ptr<ElasticMaterialBase<GridView, double>>`.
The parser always assumes the material to be heterogeneous, with the physical identity information from the grid generation process being used to map grid elements to material implementations.

The material parser accepts the following keys:
* `materials` is a list of strings that identifies subsections of the material section. Each such section describes a homogenous material.

Each material subsection accepts the following keys:
* `model` is the string that identifies the material law to be used. Currently implemented material laws are:
    - `linear` for linear isotropic elasticity theory
    - `stvenantkirchhoff` for the Saint-Venant-Kirchhoff model
    - `neohookean` for a nearly incompressible Neo-Hookean material law
* `group` is the integer that identifies the physical entity group from the GMSH file. It defaults to `0` and needs to match the `physical` parameter in the grid section.

Depending on the material law, the material subsection may contain additional keys.
The name of these keys depends on the material law as well.
For the above implemented, these are `first_lame` and `second_lame`.
For convenience, conversion from other more frequently used parameter pairs are implemented in `dune/structures/material.hh`.
E.g. the Lame parameters can be expressed through the keys `youngs_modulus` and `poisson_ratio`.

An example section for a heterogeneous material with two material is given here:

```
[material]
materials = cyto, nucleus

[material.cyto]
group = 0
youngs_modulus = 1000
poisson_ratio = 0.4
model = linear

[material.nucleus]
group = 1
youngs_modulus = 10000
poisson_ratio = 0.4
model = linear
```

TODO: Write about prestress!

### One-to-one Checker

This step is registered as `onetoone`.
It does not accept additional configuration keys.
It implements a check of the displacement mapping being one-to-one.
Add this step to your solver if you suspect your solutions to be unphysically displaced.
If this step prints failure, the calculated solution cannot be trusted.

## Utility building blocks

### Parameter definition

This step is registered as `parameter`.
It introduces a user-defined parameter into the parameter system.
The parameter will be accessible both from your custom solver steps, as well as from all your function expressions.
This way, you can specify physical parameters exactly once - removing a potential source of bugs.

The parameter step accepts the following keys:
* `name` is the string that defines the parameter name
* `datatype` is currently one of `"double"`, `"int"` or `"string"`.
* `value` is the value that this parameter is initialized with

The following example illustrates the use of a user-defined parameter in interpolation:

```
[solver]
steps = parameter, interpolation

[parameter]
name = xdispl
value = 0.5
datatype = double

[interpolation]
functions = xdispl, 0.0, 0.0
```

### Grid function probe

This step is registered as `probe`.
It adds a `GridFunctionProbe` to the solver, which evaluates the solution at a given global coordinate and prints out the result.
This is helpful when evaluating physical quantities of interest.

The probe steps accepts the following keys:
* `name` is a string that identifies the probe in the output (there might be several)
* `position` is a global coordinate given as a space-separated list of coordinates. It determines the position where the solution is supposed to be evaluated.
* `displace` is a boolean value that defaults to false. If set to true, the probe will return $`\mathbf{u}(\mathbf{x}) + \mathbf{x}`$ instead of $`\mathbf{u}(\mathbf{x})`$, which might be more insightful if the solution is a displacement field.

### Solution Transformation

This step is registered as `transformation`.
It applies a transformation to the solution.
It allows the same kinds of expression as the interpolation step, except that the additional symbols `ux`, `uy` and `uz` are available.

The transformation step accepts the following keys:
* `functions` is the list of expressions that define the transformation - one for each leaf of the function space tree.

The following block shifts a body in $`x`$-direction by adding $`1`$ to the displacement in $`x`$-direction.

```
[transformation]
functions = ux + 1.0, uy, uz
```

## Transition solver building blocks

When solving non-linear problems, it is often useful to not solve the full nonlinear problem immediately, but instead provide a sequence of problems of increasing complexity.
This section describes how such transition solvers can be implemented in dune-structures.
We distinguish variations of continuous parameter, discrete parameters and (discrete) material variations.

### Continuous parameter variation

This step is registered as `continuousvariation`.
It defines the variation of a continuous parameter, that the step introduces into the parameter system.
It always assumes the data type to be `double`.
As all the variation steps, this step holds a number of substeps to perform for each value of the varied parameter.

The continuous variation solver step accepts the following keys:
* `steps` is a list of substeps that will be performed for each value of the varied parameter
* `name` is the name of the parameter in the parameter system
* `iterations` is the number of steps to take
* `start` is the starting value of the variation parameter
* `end` is the end value of the variation parameter

The following example implements an incremental load parameter variation, where the solver step constructed for `"mysolver"` is responsible of retrieving the value of `"load"` from the parameter system (or using it through a parsed parameter function):

```
[continuousvariation]
steps = mysolver
name = load
iterations = 10
start = 0.0
end = 1000.0
```

### Discrete parameter variation

This step is registered as `discretevariation`.
It defines the variation of a discrete parameter, that the step introduces into the parameter system.
As all the variation steps, this step holds a number of substeps to perform for each value of the varied parameter.

The discrete variation solver step accepts the following keys:
* `steps` is a list of substeps that will be performed for each value of the varied parameter
* `name`is the name of the parameter in the parameter system
* `datatype` is currently one of `"double"`, `"int"` or `"string"`.
* `values` is a space-separated list of values that is to be used for this parameter

### Material parameter variation

This step is registered as `discretematerialvariation`.
It defines the variation of a discrete material parameter, that the step introduces into the parameter system.
As all the variation steps, this step holds a number of substeps to perform for each value of the varied parameter.

The discrete variation solver step accepts the following keys:
* `steps` is a list of substeps that will be performed for each value of the varied parameter
* `name`is the name of the parameter in the parameter system
* `key` is the name of the key in the material section that the parameter corresponds to
* `datatype` is currently one of `"double"`, `"int"` or `"string"`.
* `values` is a space-separated list of values that is to be used for this parameter

The following example applies a solvet first to a problem with a linear material law and then with a nonlinear one:

```
[material.mymaterial]
model = ...
...

[discretevariation]
steps = mysolver
name = modelvar
key = mymaterial.model
datatype = string
values = linear, neohookean
```

## Visualization building blocks

Visualization is integrated into the solver process.
In order to trigger the writing of visualization data, add a visualization step from the next section to your solver.

### Visualization step

This step is registered as `visualization`.
This is the main solver step whose `apply` method writes a visualization file to the disk.
For instationary problems, this step can be nested inside the time stepping solver step in order to produce an output per time step.
In order to control the datasets contained in the visualization, add substeps to this step.
Such substeps need to inherit from the base class `VisualizationStepBase<Vector>`.

The visualization step accepts the following keys:
* `instationary` is a bool specifying whether time sequences need to be written
* `name` is the basename of the output file(s)
* `path` is a relative directory where to place the output. If the directory does not exist, it is created.

### Solution visualization step

This step is registered as `vis_solution`.
It needs to be added as a child step to a visualization node.
It visualizes the current solution and accepts no configuration values.

### MPI rank visualization step

This step is registered as `vis_mpirank`.
It needs to be added as a child step to a visualization node.
It visualizes the MPI ranks in a distributed mesh and accepts no configuration values.

### Von-Mises Stress visualization step

This step is registered as `vis_vonmises`.
It needs to be added as a child step to a visualization node.
It visualizes the Von-Mises stress of a given displacement field and accepts no configuration values.

### Physical entity data visualization step

This step is registered as `vis_physicalentity`.
It needs to be added as a child step to a visualization node.
It visualizes the physical entity information from the GMSH file and accepts no configuration values.

### Fibre distance visualization step

This step is registered as `vis_fibredistance`.
It needs to be added as a child step to a visualization node.
It visualizes the distance to a fibre as used in fibre-aligned prestressing.
This is much more a debugging tool than something you would usually want to use.

It accepts the following keys:
* `key` is the key of the fibre configuration section from the grid section

# How to...

here are some guidelines for likely extensions of the dune-structures code base.

## ... add a material law

- Have a look at the file `python/dune/structures/material.py`.
- Add a class to that file that inherits from `MaterialLawBase` and implements its fields
- Add the ID and material names you added in python to the `std::map` data structures in `dune/structures/material.hh`
- Register your material for code generation by adding it to the list of materials in `operators/operators.ufl`

## ... add a solver step

* Write a class that inherits from `TransitionSolverStepBase<Vector>` from the file `dune/structures/solversteps/base.hh`.
* Decide whether a more specialized base class is useful for you:
    - `StepCollectionStep<Vector>` implements a step that accepts substeps
    - `WrapperStep<Vector>` implements a step that wraps another step
    - `VisualizationStepBase<Vector>` implements a step that adds data to visualization
* If you want to construct your step from a configuration file, you must register a construction method like seen in `dune/structures/solverconstruction.hh`. This can either be done in that file or you can call the `registerStep` method of the `ConstructionContext` object directly or by doing this in the constructor of the construction object.
* If your solver step requires parameter functions, have a look at the `get_callable` and `get_callable_array` functions from `dune/structures/solverconstruction.hh`.

## ... do XYZ

Write an email to [dominic.kempf@iwr.uni-heidelberg.de](mailto:dominic.kempf@iwr.uni-heidelberg.de).