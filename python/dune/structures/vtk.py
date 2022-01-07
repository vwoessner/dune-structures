import os
from collections.abc import Mapping

import numpy as np
from vtk import vtkXMLGenericDataObjectReader, vtkUnstructuredGrid


class VTKUnstructuredGridReader:
    _grid = None  # The grid accessed by this reader

    def __init__(self, filepath):
        # Check file path
        filepath = os.path.abspath(filepath)
        if not os.path.isfile(filepath):
            raise FileNotFoundError("File not found: {}".format(filepath))

        # Set up the reader
        reader = vtkXMLGenericDataObjectReader()
        reader.SetFileName(filepath)
        reader.Update()
        self._grid = vtkUnstructuredGrid.SafeDownCast(reader.GetOutput())


class VTKVertexReader(VTKUnstructuredGridReader, Mapping):
    def __init__(self, filepath, ensure_anticlockwise_points=False):
        # Load the grid
        super().__init__(filepath)

        # Retrieve the points and determine unique ones
        self._points = np.array(
            [self._grid.GetPoint(idx) for idx in range(self._grid.GetNumberOfPoints())]
        )
        self._points, self._unique_idx, self._unique_idx_inv = np.unique(
            self._points, axis=0, return_index=True, return_inverse=True
        )

        # Get connectivity information
        cell_array = self._grid.GetCells()
        connect_array = cell_array.GetConnectivityArray()
        self._connectivity = np.array(
            [
                connect_array.GetTuple(idx)
                for idx in range(cell_array.GetNumberOfConnectivityIds())
            ],
            dtype=int,
        )[self._unique_idx_inv]

        # Print offsets
        # offsets = cell_array.GetOffsetsArray()
        # print([offsets.GetTuple(idx) for idx in range(cell_array.GetNumberOfOffsets())])

        # NOTE: Assuming triangles!
        self._connectivity = np.reshape(self._connectivity, (-1, 3))

        # Order connectivity data
        if ensure_anticlockwise_points:
            for idx in range(self._connectivity.shape[0]):
                idx_sort = np.argsort(
                    self._center_arc(self._points[self._connectivity[idx]])
                )
                self._connectivity[idx] = self._connectivity[idx, idx_sort]

        # Prepare access to data
        self._point_data = self._grid.GetPointData()

    def _center_arc(self, coordinates):
        center = np.mean(coordinates, axis=0)
        rad_vec = coordinates - center
        return np.arctan2(rad_vec[..., 1], rad_vec[..., 0])

    def __len__(self):
        return self._point_data.GetNumberOfArrays()

    def __iter__(self):
        for idx in range(self.__len__()):
            yield self._point_data.GetArrayName(idx)

    def __getitem__(self, key):
        if self._point_data.HasArray(key) == 1:
            data_array = self._point_data.GetArray(key)
            return np.squeeze(
                np.array(
                    [
                        data_array.GetTuple(idx)
                        for idx in range(self._grid.GetNumberOfPoints())
                    ]
                )
            )[self._unique_idx]
        else:
            raise KeyError("No dataset found with name: {}".format(key))

    @property
    def points(self):
        return self._points

    @property
    def connectivity(self):
        return self._connectivity

    @property
    def unique_idx(self):
        return self._unique_idx

    @property
    def unique_idx_inv(self):
        return self._unique_idx_inv


class VTKCellReader(VTKUnstructuredGridReader, Mapping):
    def __init__(self, filepath):
        # Load the grid
        super().__init__(filepath)

        # Retrieve the cell positions
        self._cell_points = []
        self._cell_centers = []
        for idx in range(self._grid.GetNumberOfCells()):
            cell = self._grid.GetCell(idx)
            points = [
                cell.GetPoints().GetPoint(i) for i in range(cell.GetNumberOfPoints())
            ]
            self._cell_points.append(points)
            self._cell_centers.append(np.mean(points, axis=0))

        # Prepare access to data
        self._cell_data = self._grid.GetCellData()

    def __len__(self):
        return self._cell_data.GetNumberOfArrays()

    def __iter__(self):
        for idx in range(self.__len__()):
            yield self._cell_data.GetArrayName(idx)

    def __getitem__(self, key):
        if self._cell_data.HasArray(key) == 1:
            data_array = self._cell_data.GetArray(key)
            return np.squeeze(
                np.array(
                    [
                        data_array.GetTuple(idx)
                        for idx in range(self._grid.GetNumberOfCells())
                    ]
                )
            )
        else:
            raise KeyError("No dataset found with name: {}".format(key))

    @property
    def cell_points(self):
        """Return the coordinates of the points spanning the cells. The
        dimension of this array is 3: Index of the cell, index of the point, three
        coordinates. The object is returned as list that can be directly inserted into a
        ``matplotlib.poly.PolyCollection``."""
        return self._cell_points

    @property
    def cell_centers(self):
        return np.asarray(self._cell_centers)
