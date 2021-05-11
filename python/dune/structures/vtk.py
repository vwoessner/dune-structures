import os
from collections.abc import Mapping

import numpy as np
from vtk import vtkXMLGenericDataObjectReader, vtkUnstructuredGrid


class VTKVertexReader(Mapping):
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
            ]
        )[self._unique_idx_inv]

        # NOTE: Assuming triangles!
        self._connectivity = np.reshape(self._connectivity, (-1, 3))

        # Prepare access to data
        self._point_data = self._grid.GetPointData()

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
