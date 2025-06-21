import math

import numpy as np
import pyvista
from scipy.spatial.transform import Rotation
import vtk
from vtk.util import numpy_support # type: ignore
from vtkbool import vtkBool

from Trace.AbstractTrace import AbstractTrace
from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
import Geometry.geometry_tools as geo
from Trace.PipeShape import PipeShape
from Geometry.LineSegment import LineSegment
import tool.vtk_tools as vt

class VtkPointGroup(PntInc):
    def __init__(self, xyz_points: np.ndarray, vtk_indices: list[int] = None):
        if vtk_indices is None:
            vtk_indices = [None for i in range(xyz_points.shape[0])]

        self.xyz_points = xyz_points
        self.vtk_indices = vtk_indices

    @property
    def vtk_idx_0(self) -> int | None:
        return self.vtk_indices[0]

    @property
    def vtk_idx_n(self) -> int | None:
        return self.vtk_indices[-1]

    def add_missing_vtk_points(self, polydata: vtk.vtkPolyData):
        vtk_points: vtk.vtkPoints = polydata.GetPoints()

        for i in range(len(self.xyz_points)):
            if self.vtk_indices[i] is None:
                xyz_points: np.ndarray = self.xyz_points[i]
                vtk_points.InsertNextPoint(xyz_points.tolist())
                self.vtk_indices[i] = vtk_points.GetNumberOfPoints()-1

    def inc_vtk_indicies(self, start: int, cnt: int):
        for i in range(len(self.vtk_indices)):
            if self.vtk_indices[i] >= start:
                self.vtk_indices[i] += cnt
    
    def __len__(self):
        return self.xyz_points.shape[0]