import pytest

import math
import numpy as np

from dune.structures import VTKVertexReader

FILE = "test.vtu"

# --- Fixtures ---


@pytest.fixture
def vtk_reader():
    """Create a VTKReader for test file"""
    return VTKVertexReader(FILE)


# --- Tests ---


def test_init(vtk_reader):
    """Test initialization of the reader"""
    assert len(vtk_reader) == 2
    assert "vectors" in vtk_reader
    assert "scalars" in vtk_reader

    with pytest.raises(KeyError):
        vtk_reader["does-not-exist"]

    with pytest.raises(FileNotFoundError):
        VTKVertexReader("does-not-exist.vtu")


def test_points(vtk_reader):
    """Check if unique points are correctly detected"""
    assert len(vtk_reader.points) == 4
    assert vtk_reader.connectivity.shape == (2, 3)
    assert len(vtk_reader.unique_idx) == 4
    assert len(vtk_reader.unique_idx_inv) == 6


def test_data(vtk_reader):
    """Verify data"""
    scalars = vtk_reader["scalars"]
    assert scalars.shape == (4,)
    assert math.isclose(float(np.mean(scalars)), 0.5)

    vectors = vtk_reader["vectors"]
    assert vectors.shape == (4, 3)
    assert math.isclose(float(np.mean(vectors)), 0.5)
