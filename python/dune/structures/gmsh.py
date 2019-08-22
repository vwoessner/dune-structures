""" Create idealized cell geometries for meshing with GMSH - through a convenient Python interface """

import hashlib
import meshio
import os.path
import pygmsh
import sys
import yaml
import numpy as np
import shutil
import subprocess
import pytools as pt


def _parse_vec(x):
    if isinstance(x, str):
        x = x.split()
    return np.array(list(float(i) for i in x))


class GMSHError(Exception):
    pass


@pt.memoize
def get_gmsh_version(gmshexec):
    out = subprocess.check_output([gmshexec, "--version"], stderr=subprocess.STDOUT).strip().decode("utf8")

    def mysplit(s, delims):
        if len(delims) == 0:
            yield s
            return

        for part in s.split(delims[0]):
            yield from mysplit(part, delims[1:])

    return mysplit(out, [".", "-"])


@pt.memoize
def get_gmsh_build_flags(gmshexec):
    out = subprocess.check_output([gmshexec, "-info"], stderr=subprocess.STDOUT).strip().decode("utf8")
    out = out.split("\n")
    for line in out:
        if line.startswith("Build options"):
            line = line[line.find(":") + 1:]
            return [o.strip() for o in line.split()]

    raise GMSHError("The output format of 'gmsh -info' was unexpected!")


def check_gmsh_build_flag(gmshexec, flag):
    return flag in get_gmsh_build_flags(gmshexec)


def generate_cell_mesh(config, mshfile, gmshexec="gmsh"):
    """ The entry point for the creation of a cell mesh """
    geo = pygmsh.opencascade.Geometry()

    geo.add_comment("Generated for gmsh version {}".format(".".join(get_gmsh_version(gmshexec))))

    # A dictionary with material information
    material_to_geo = {}

    def add_shape(config):
        shape = config.get("shape", "box")
        if shape == "box":
            size = _parse_vec(config.get("size", [1.0, 1.0, 1.0]))
            lowerleft = _parse_vec(config.get("lowerleft", -0.5 * size))
            return geo.add_box(lowerleft, size)
        elif shape == "sphere":
            radius = config.get("radius", 1.0)
            center = _parse_vec(config.get("center", [0.0, 0.0, 0.0]))
            return geo.add_ball(center, radius)
        elif shape == "round":
            radius = config.get("radius", 1.0)
            center = _parse_vec(config.get("center", [0.0, 0.0, 0.6]))
            cutoff = config.get("cutoff", None)
            return geo.add_ball(center, radius, x0=cutoff)
        elif shape == "spread":
            radius = config.get("radius", 1.0)
            height = config.get("height", 1.0)
            # The slope parameter is from (0,1) and describes the position of the third control point of the BSpline
            slope = config.get("slope", 0.5)

            # Necessary points
            b0 = geo.add_point([0.0, 0.0, height])
            b1 = geo.add_point([0.5 * radius, 0.0, height])
            b2 = geo.add_point([radius, 0.0, 0.0])
            b3 = geo.add_point([0.0, 0.0, 0.0])

            # Construct a planar surface that can be rotated
            spline = geo.add_bspline([b0, b1, b2])
            loop = geo.add_line_loop([spline, geo.add_line(b2, b3), geo.add_line(b3, b0)])
            surf = geo.add_plane_surface(loop)

            cell = geo.extrude(surf,
                               rotation_axis=[0.0, 0.0, 1.0],
                               point_on_axis=[0.0, 0.0, 0.0],
                               angle=2.0*np.pi
                               )

            # Only return the volume ID
            return cell[1]
        elif shape == "ellipsoid":
            center = _parse_vec(config.get("center", [0.0, 0.0, 0.0]))
            radii = _parse_vec(config.get("radii", [2.0, 1.0, 1.0]))
            cutoff = config.get("cutoff", None)
            return geo.add_ellipsoid(center, radii, x0=cutoff)
        else:
            raise NotImplementedError

    def add_material(geoobj, config):
        physical = config.get("physical", None)
        if physical is not None:
            material_to_geo.setdefault(physical, [])
            material_to_geo[physical].append(geoobj)

    # Add the cytoplasma
    cytoconfig = config.get("cytoplasm", {})
    cyto = add_shape(cytoconfig)

    # Maybe add a nucleus
    nucleusconfig = config.get("nucleus", {"enabled": False})
    if nucleusconfig.get("enabled", True):
        # Add the nucleus shape and intersect it with the cytoplasma
        nucleus = add_shape(nucleusconfig)
        nucleus = geo.boolean_intersection([cyto, nucleus], delete_first=False, delete_other=True)

        add_material(nucleus, nucleusconfig)

    # Maybe add some fibres
    fibres = []
    fibreconfig = config.get("fibres", [])
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
        add_material(fibre, fconfig)

        # Store the fibre object for later
        fibres.append(fibre)

    # Implement mesh widths
    geo.add_raw_code("Characteristic Length{{ PointsOf{{ Volume{{:}}; }} }} = {};".format(cytoconfig.get("meshwidth", 0.1)))
    if nucleusconfig.get("enabled", True):
        meshwidth = nucleusconfig.get("meshwidth", cytoconfig.get("meshwidth", 0.1))
        geo.add_raw_code("Characteristic Length{{ PointsOf{{ Volume{{{}}}; }} }} = {};".format(nucleus.id, meshwidth))
    for i, fibre in enumerate(fibres):
        geo.add_raw_code("Characteristic Length{{ PointsOf{{ Volume{{{}}}; }} }} = {};".format(fibre.id, fibreconfig[i].get("meshwidth", 0.02)))

    # The cytoplasma is defined by the outer shape minus all inclusions
    cytogeos = [cyto]
    if nucleusconfig.get("enabled", True):
        cytogeos.append(nucleus)
    cyto = geo.boolean_fragments(cytogeos, fibres, delete_other=False)
    add_material(cyto, cytoconfig)

    # Add the collected physical group information
    for physical, geos in material_to_geo.items():
        geo.add_physical(geos, physical)

    # The default meshing algorithm creates spurious elements on the sphere surface
    if cytoconfig.get("shape", "box") != "box":
        geo.add_raw_code("Mesh.Algorithm = 6;")

    # Maybe skip mesh generation if we have a cache hit
    cachefile = "{}.msh".format(hashlib.sha256("\n".join(geo.get_code()).encode()).hexdigest())
    if os.path.exists(cachefile):
        print("A cached version of this mesh was found!")
        shutil.copyfile(cachefile, mshfile)
        return

    # Maybe output a geofile. We do so before calling gmsh to use it for debugging
    exportconfig = config.get("export", {})
    geoconfig = exportconfig.get("geo", {})
    if geoconfig.get("enabled", False):
        filename = geoconfig.get("filename", mshfile)
        filename = "{}.geo".format(os.path.splitext(mshfile)[0])
        with open(filename, 'w') as f:
            f.write(geo.get_code() + "\n")

    # Finalize the grid generation
    try:
        mesh = pygmsh.generate_mesh(geo, gmsh_path=gmshexec)
    except AssertionError:
        raise GMSHError("Gmsh failed. Check if nucleus intersects fibres (which is not supported)!")

    # Export this mesh into several formats as requested
    mshconfig = exportconfig.get("msh", {})
    meshio.write(mshfile, mesh, file_format="gmsh2-ascii")
    shutil.copyfile(mshfile, cachefile)

    vtkconfig = exportconfig.get("vtk", {})
    if vtkconfig.get("enabled", False):
        filename = vtkconfig.get("filename", mshfile)
        filename = "{}.vtk".format(os.path.splitext(mshfile)[0])
        # This throws a warning, but the physical groups only work in ASCII mode
        meshio.write(filename, mesh, write_binary=False)


def entrypoint_generate_mesh():
    # Read the command line arguments to this script
    config = yaml.safe_load(open(sys.argv[1], 'r'))
    mshfile = sys.argv[2]
    gmshexec = "gmsh"
    if len(sys.argv) > 3:
        gmshexec = sys.argv[3]

    generate_cell_mesh(config, mshfile, gmshexec)
