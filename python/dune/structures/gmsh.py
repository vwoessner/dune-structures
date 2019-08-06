""" Create idealized cell geometries for meshing with GMSH - through a convenient Python interface """

import meshio
import pygmsh
import sys
import yaml


def _parse_vec(x):
    if isinstance(x, str):
        x = x.split()
    return list(float(i) for i in x)


def generate_cell_mesh(config, mshfile):
    """ The entry point for the creation of a cell mesh """
    geo = pygmsh.opencascade.Geometry()

    # Add the cytoplasma
    cytoconfig = config.get("cytoplasm", {})
    shape = cytoconfig.get("shape", "box")
    if shape == "box":
        size = _parse_vec(cytoconfig.get("size", [1.0, 1.0, 1.0]))
        lowerleft = list(-0.5 * i for i in size)
        cyto = geo.add_box(lowerleft, size)
    elif shape == "sphere":
        raise NotImplementedError
    elif shape == "round":
        raise NotImplementedError
    else:
        raise NotImplementedError

    # Maybe add a nucleus
    nucleusconfig = config.get("nucleus", {})
    if nucleusconfig.get("enabled", False):
        raise NotImplementedError

    # Maybe add some fibres
    fibres = []
    fibreconfig = config.get("fibres", [])
    for fconfig in fibreconfig:
        # Parse fiber data
        offset = _parse_vec(fconfig.get("offset", [-1.0, 0.0, 0.0]))
        orientation = _parse_vec(fconfig.get("orientation", [2.0, 0.0, 0.0]))
        radius = float(fconfig.get("radius", 0.1))
        
        # Add a raw cylinder
        fibre = geo.add_cylinder(offset, orientation, radius)

        # And intersect it with the cytoplasma
        fibres.append(geo.boolean_intersection([cyto, fibre], delete_first=False, delete_other=True))

    # Post process fibre generation
    geo.boolean_fragments([cyto], fibres, delete_other=False)

    # Implement mesh widths
    geo.add_raw_code("Characteristic Length{{ PointsOf{{ Volume{{:}}; }} }} = {};".format(cytoconfig.get("meshwidth", 0.1)))
    for i, fibre in enumerate(fibres):
        geo.add_raw_code("Characteristic Length{{ PointsOf{{ Volume{{{}}}; }} }} = {};".format(fibre.id, fibreconfig[i].get("meshwidth", 0.02)))
    
    print(geo.get_code())

    # Finalize the grid generation
    mesh = pygmsh.generate_mesh(geo)
    meshio.write(mshfile, mesh, write_binary=False)
    meshio.write(mshfile + ".vtk", mesh)


def entrypoint_generate_mesh():
    # Read the command line arguments to this script
    config = yaml.safe_load(open(sys.argv[1], 'r'))
    mshfile = sys.argv[2]

    generate_cell_mesh(config, mshfile)
