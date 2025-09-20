import numpy as np
import pyvista
from scipy.spatial.transform import Rotation
import vtk
from vtk.util import numpy_support # type: ignore
#from vtkbool import vtkBool

from Trace.AbstractTrace import AbstractTrace
from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
import Geometry.geometry_tools as geo
from Geometry.PathCorner import PathCorner
from Geometry.LineSegment import LineSegment
from Trace.VtkPointGroup import VtkPointGroup
import tool.vtk_tools as vt

class TraceCorner(PathCorner, PntInc):
    """
    A class representing a corner in the trace where two line segments meet.

    Parameters
    ----------
    parent : AbstractTrace
        The parent instance that contains this instance.
    segments : tuple[LineSegment, LineSegment]
        The two segments that this corner is found on. They should share one of the end points.
    bend_radius : float
        How large of a radius the bend should be. The actual radius used will be the maximum of this radius and the trace shape radius.

    Attributes
    ----------
    parent : AbstractTrace
        The parent instance that contains this instance.
    xypntidx_to_vtkidx : list[VtkPointGroup]
        The VTK indicies per xy point index, indicating which VTK points correspond to each centerline point in the corner.
        """
    def __init__(self, parent: AbstractTrace, segments: tuple[LineSegment, LineSegment], bend_radius: float):
        self.parent: AbstractTrace = None  # Type hint override
        super().__init__(parent, segments, bend_radius)

        self.xypntidx_to_vtkidx: list[VtkPointGroup] = []
        """ The VTK indicies per xy point. """

    @property
    def bend_radius(self) -> float:
        """
        Overrides the PathCorner bend_radius to ensure it is at least as large as the shape's radius.
        """
        # Get the diameter from the trace's shape and convert to radius for comparison.
        shape_radius = self.parent.shape.diameter / 2
        return max(super().bend_radius, shape_radius)

    @property
    def n_segments(self) -> int:
        """ The number of segments in this corner. """
        return len(self.get_center_line_points()) - 1

    @property
    def n_points(self) -> int:
        """ The total number of center line points for this corner. """
        return len(self.get_center_line_points())

    @property
    def vtk_idx_0(self) -> int:
        """ Index in the vtk points where the first loop of the trace corner starts. """
        return self.xypntidx_to_vtkidx[0].vtk_idx_0

    def get_vtk_group(self, xy_idx: int) -> VtkPointGroup:
        """
        Gets the list of VTK indicies that make up the loop at
        the given center line xy point index or coordinate.

        Parameters
        ----------
        xy_idx : int
            The index for a specific center line point in this corner.

        Returns
        -------
        VtkPointGroup
            The VTK indicies corresponding to the given center line point.
        """
        # build up the xy points list
        while len(self.xypntidx_to_vtkidx) < xy_idx + 1:
            xy_point_idx = len(self.xypntidx_to_vtkidx)
            (x, y), angle = self.get_center_line_points()[xy_point_idx]
            loop_xyz_points = self.parent.shape.oriented_points(angle, (x, y))
            vtk_group = VtkPointGroup(np.array(loop_xyz_points))
            self.xypntidx_to_vtkidx.append(vtk_group)

        return self.xypntidx_to_vtkidx[xy_idx]

    def to_vtk(self, polydata: vtk.vtkPolyData):
        """ Adds this corner into a VTK PolyData structure. """
        # add any missing vtk points, for each xy centerline point
        for xy_idx in range(self.n_points):
            vtk_group = self.get_vtk_group(xy_idx)
            vtk_group.add_missing_vtk_points(polydata)
            vt.calculate_point_normals(polydata, vtk_group.vtk_idx_0, vtk_group.vtk_idx_n)

        # connect the adjacent loops with quads
        for a in range(self.n_points - 1):
            va = self.get_vtk_group(a)
            vb = self.get_vtk_group(a + 1)
            vt.adjoin_with_quads(polydata, va.vtk_idx_0, vb.vtk_idx_0, len(va))

    def inc_vtk_indicies(self, start: int, cnt: int):
        for xy_pnt_idx in range(len(self.xypntidx_to_vtkidx)):
            vtk_indicies = self.xypntidx_to_vtkidx[xy_pnt_idx]
            for vtkidx_idx in range(len(vtk_indicies)):
                if vtk_indicies[vtkidx_idx] >= start:
                    vtk_indicies[vtkidx_idx] += cnt