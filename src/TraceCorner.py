import numpy as np
import vtk

from AbstractTrace import AbstractTrace
from AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
import geometry as geo
from PathCorner import PathCorner
from Segment import Segment

class TraceCorner(PathCorner, PntInc):
    def __init__(self, parent: AbstractTrace, segments: tuple[Segment, Segment], bend_radius: float):
        """
        Parameters
        ----------
        parent : AbstractTrace
            The parent instance that contains this instance.
        segments : tuple[Segment, Segment]
            The two segments that this corner is found on. They should
            share one of the end points.
        bend_radius : float
            How large of a radius the bend should be. The actual radius
            used will be the maximum of this radius and the trace
            shape radius.
        """
        self.parent: AbstractTrace = None # type hint override
        super(PathCorner, self).__init__(parent, segments, bend_radius)

        self.loop_a_idx: int = None
        """ Index in the vtk points where the first loop of the trace corner starts. """
        self.loop_b_idx: int = None
        """ Index in the vtk points where the last loop of the trace corner starts. """

    @property
    def bend_radius(self):
        # overrides the PathCorner bend_radius
        shape_radius = self.parent.shape.diameter / 2
        return max(super().bend_radius, shape_radius)

    def inc_vtk_indicies(self, start: int, cnt: int):
        if self.loop_a_idx >= start:
            self.loop_a_idx += cnt
        if self.loop_b_idx >= start:
            self.loop_b_idx += cnt

    def to_vtk(self, polydata: vtk.vtkPolyData):
        vtk_points = polydata.GetPointData()
        vtk_cells = polydata.GetCellData()

        # get some properties of the two segments
        xy_point = self.parent.xy_points[self.segments[1]]

        # get properties of the arc between the two segments
        angle_diff, mid_angle, arc_length, arc_center = self.get_arc_properties()

        # 