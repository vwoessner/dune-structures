""" Create idealized cell geometries for meshing with GMSH - through a convenient Python interface """

import meshio
import os.path
import pygmsh
import sys
import yaml


def _parse_vec(x):
    if isinstance(x, str):
        x = x.split()
    return list(float(i) for i in x)


def add_shape(geo, config):
    shape = config.get("shape", "box")
    if shape == "box":
        size = _parse_vec(config.get("size", [1.0, 1.0, 1.0]))
        lowerleft = list(-0.5 * i for i in size)
        return geo.add_box(lowerleft, size)
    elif shape == "sphere":
        radius = config.get("radius", 1.0)
        return geo.add_ball([0.0, 0.0, 0.0], radius)
    elif shape == "round":
        raise NotImplementedError
    else:
        raise NotImplementedError


def generate_cell_mesh(config, mshfile):
    """ The entry point for the creation of a cell mesh """
    geo = pygmsh.opencascade.Geometry()

    # Add the cytoplasma
    cytoconfig = config.get("cytoplasm", {})
    cyto = add_shape(geo, cytoconfig)

    # Maybe add a nucleus
    nucleusconfig = config.get("nucleus", {"enabled": False})
    if nucleusconfig.get("enabled", True):
        # Add the nucleus shape
        nucleus = add_shape(geo, nucleusconfig)

        #TODO Physical data
        nucleus = geo.boolean_intersection([cyto, nucleus], delete_first=False, delete_other=True)

    # Maybe add some fibres
    fibres = []
    fibreconfig = config.get("fibres", [])
    physical_to_fibre = {}
    for i, fconfig in enumerate(fibreconfig):
        # Parse fiber data
        offset = _parse_vec(fconfig.get("offset", [-1.0, 0.0, 0.0]))
        orientation = _parse_vec(fconfig.get("orientation", [2.0, 0.0, 0.0]))
        radius = float(fconfig.get("radius", 0.1))
        
        # Add a raw cylinder
        fibre = geo.add_cylinder(offset, orientation, radius)

        # And intersect it with the cytoplasma
        fibre = geo.boolean_intersection([cyto, fibre], delete_first=False, delete_other=True)

        # Add physical information to this fibre
        fibre_physical = fconfig.get("physical", None)
        if fibre_physical is not None:
            physical_to_fibre.setdefault(fibre_physical, [])
            physical_to_fibre[fibre_physical].append(fibre)

        # Store the fibre object for later
        fibres.append(fibre)

    # Post process fibre generation
    frags = geo.boolean_fragments([cyto, nucleus], fibres, delete_other=False)

    # Implement mesh widths
    geo.add_raw_code("Characteristic Length{{ PointsOf{{ Volume{{:}}; }} }} = {};".format(cytoconfig.get("meshwidth", 0.1)))
    for i, fibre in enumerate(fibres):
        geo.add_raw_code("Characteristic Length{{ PointsOf{{ Volume{{{}}}; }} }} = {};".format(fibre.id, fibreconfig[i].get("meshwidth", 0.02)))

    # Maybe add physical group information
    for physical, fibres in physical_to_fibre.items():
        geo.add_physical(fibres, physical)

    nucleus_physical = nucleusconfig.get("physical", None)
    if nucleusconfig.get("enabled", True) and nucleus_physical is not None:
        geo.add_physical([nucleus], nucleus_physical)

    cyto_physical = cytoconfig.get('physical', None)
    if cyto_physical is not None:
        geo.add_physical([frags], cyto_physical)

    # The default meshing algorithm creates spurious elements on the sphere surface
    if cytoconfig.get("shape", "box") != "box":
        geo.add_raw_code("Mesh.Algorithm = 5;")

    # Finalize the grid generation
    mesh = pygmsh.generate_mesh(geo)

    # Export this mesh into several formats as requested
    exportconfig = config.get("export", {})

    mshconfig = exportconfig.get("msh", {})
    meshio.write(mshfile, mesh, write_binary=False)

    vtkconfig = exportconfig.get("vtk", {})
    if vtkconfig.get("enabled", False):
        filename = vtkconfig.get("filename", mshfile)
        filename = "{}.vtk".format(os.path.splitext(mshfile)[0])
        # This throws a warning, but the physical groups only work in ASCII mode
        meshio.write(filename, mesh, write_binary=False)

    geoconfig = exportconfig.get("geo", {})
    if geoconfig.get("enabled", False):
        filename = geoconfig.get("filename", mshfile)
        filename = "{}.geo".format(os.path.splitext(mshfile)[0])
        with open(filename, 'w') as f:
            f.write(geo.get_code() + "\n")


def entrypoint_generate_mesh():
    # Read the command line arguments to this script
    config = yaml.safe_load(open(sys.argv[1], 'r'))
    mshfile = sys.argv[2]

    generate_cell_mesh(config, mshfile)
