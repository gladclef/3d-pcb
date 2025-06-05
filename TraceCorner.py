import numpy as np
import vtk

from AbstractTrace import AbstractTrace
from AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc

class TraceCorner(PntInc):
    def __init__(self, parent: AbstractTrace, segments: tuple[int, int, int], bend_radius: float):
        """
        Parameters
        ----------
        parent : AbstractTrace
            The parent instance that contains this instance.
        segments : tuple[int, int, int]
            The two segments that this corner is found on, as the three
            xy point ids in the parent instance.
        bend_radius : float
            How large of a radius wire bends should be. The maximum of
            this value and the trace's pipe shape will be used as the
            actual bend radius.
        """
        self.parent = parent
        self.segments = segments
        self._bend_radius = bend_radius

        self.loop_a_idx: int = None
        """ Index in the vtk points where the first loop of the trace corner starts. """
        self.loop_b_idx: int = None
        """ Index in the vtk points where the last loop of the trace corner starts. """

    @property
    def bend_radius(self):
        shape_radius = self.parent.shape.diameter / 2
        return max(self._bend_radius, shape_radius)
    
    @bend_radius.setter
    def bend_radius(self, val: float):
        self._bend_radius = val

    def _get_arc_properties(self, angle_a: float, angle_b: float, xy_point: tuple[float, float]) -> tuple[float, float, tuple[float, float]]:
        # compute the arc angles
        angle_diff = 0 # the angle to add to angle_a to get to angle_b
        _mid_angle = 0 # the angle of the line that is halfway between angle_a and angle_b
        angle_diff = np.abs(angle_b - angle_a)
        angle_diff = angle_diff if angle_diff > np.pi else 2*np.pi - angle_diff
        mid_angle = angle_diff
        mid_rise_over_run = np.sin(mid_angle) / np.cos(mid_angle)
        arc_center_angle = -1/mid_rise_over_run # the angle

        # compute the length of the arc necessary based on the given bend radius and the two angles a and b
        arc_length = (angle_diff / (2 * np.pi)) * 2 * np.pi * self.bend_radius

        # determine the location of the center of the arc such that it goes through the xy_point
        mid_angle = angle_a + (angle_a + angle_b) / 2
        center_x = xy_point[0] + self.bend_radius * np.cos(mid_angle)
        center_y = xy_point[1] + self.bend_radius * np.sin(mid_angle)

        expected_angle_diff = np.pi*2 - abs(angle_a - angle_b)
        halfcirc_angle_diff = abs(angle_a - angle_b)
        expected_arc_length = (halfcirc_angle_diff / (2 * np.pi)) * (2 * np.pi * self.instance.bend_radius)
        mid_angle = (angle_a + angle_b) / 2
        rise_over_run = np.sin(mid_angle) / np.cos(mid_angle)
        arc_center_angle = -1/rise_over_run
        expected_center_x = xy_point[0] + self.instance.bend_radius * np.cos(arc_center_angle)
        expected_center_y = xy_point[1] + self.instance.bend_radius * np.sin(arc_center_angle)

        return angle_diff, arc_length, (center_x, center_y)

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
        angle_a = self.parent.get_setment_angle((self.segments[0], self.segments[1]))
        angle_b = self.parent.get_setment_angle((self.segments[1], self.segments[2]))

        # get properties of the arc between the two segments
        angle_diff, arc_length, xy_arc_center = self._get_arc_properties(angle_a, angle_b, xy_point)

        # 