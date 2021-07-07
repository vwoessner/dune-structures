import pytest

import math
import numpy as np

from dune.structures import VTKVertexReader, VTKCellReader

FILE = "test.vtu"

# --- Fixtures ---


@pytest.fixture
def vtk_vertex_reader():
    """Create a VTKReader for test file"""
    return VTKVertexReader(FILE)


@pytest.fixture
def vtk_cell_reader():
    """Create a VTKReader for test file"""
    return VTKCellReader(FILE)


@pytest.fixture
def vtk_readers(vtk_vertex_reader, vtk_cell_reader):
    """Merge both readers into a dict"""
    return {"vertex": vtk_vertex_reader, "cell": vtk_cell_reader}


# --- Tests ---


def test_init_throw(vtk_readers):
    """Test initialization of the readers"""
    for reader in vtk_readers.values():
        with pytest.raises(KeyError):
            reader["does-not-exist"]

    with pytest.raises(FileNotFoundError):
        VTKVertexReader("does-not-exist.vtu")
    with pytest.raises(FileNotFoundError):
        VTKCellReader("does-not-exist.vtu")


# Vertex Reader


def test_vertex_init(vtk_vertex_reader):
    assert len(vtk_vertex_reader) == 2
    assert "vectors" in vtk_vertex_reader
    assert "scalars" in vtk_vertex_reader
    assert "scalars2" not in vtk_vertex_reader


def test_vertex_points(vtk_vertex_reader):
    """Check if unique points are correctly detected"""
    assert len(vtk_vertex_reader.points) == 4
    assert vtk_vertex_reader.connectivity.shape == (2, 3)
    assert len(vtk_vertex_reader.unique_idx) == 4
    assert len(vtk_vertex_reader.unique_idx_inv) == 6


def test_vertex_data(vtk_vertex_reader):
    """Verify data"""
    scalars = vtk_vertex_reader["scalars"]
    assert scalars.shape == (4,)
    assert math.isclose(float(np.mean(scalars)), 0.5)

    vectors = vtk_vertex_reader["vectors"]
    assert vectors.shape == (4, 3)
    assert math.isclose(float(np.mean(vectors)), 0.5)


# Cell Reader


def test_cell_init(vtk_cell_reader):
    assert len(vtk_cell_reader) == 1
    assert "scalars2" in vtk_cell_reader
    assert "vectors" not in vtk_cell_reader


def test_cell_points(vtk_cell_reader):
    cell_points = vtk_cell_reader.cell_points
    assert isinstance(cell_points, list)
    assert len(cell_points) == 2  # 2 cells
    cell_points = np.array(cell_points)
    assert cell_points.shape == (2, 3, 3)  # 2 cells, 3 points, 3 coords

    assert vtk_cell_reader.cell_centers.shape == (2, 3)


def test_cell_data(vtk_cell_reader):
    scalars = vtk_cell_reader["scalars2"]
    assert scalars.shape == (2,)
    assert math.isclose(float(np.mean(scalars)), 2.0)
