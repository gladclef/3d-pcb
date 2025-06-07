import numpy as np
import vtk

from AbstractTrace import AbstractTrace
from AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
import geometry as geo
from Segment import Segment

class TraceCorner(PntInc):
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
            How large of a radius wire bends should be. The maximum of
            this value and the trace's pipe shape will be used as the
            actual bend radius.
        """
        # verify that the segments share an end point, and that the shared point is in the middle
        if segments[0].xy1 == segments[1].xy0:
            if segments[1].xy1 == segments[0].xy0:
                raise ValueError("The segments must share only one common end point.")
        elif segments[1].xy1 == segments[0].xy0:
            segments = [segments[1], segments[0]]
        else:
            raise ValueError("The segments must share a common end point.")

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

    def get_arc_seg_intersections(self):
        """
        Find the location where the arc intersects both line segments.

        To do this we:
            1. assume both lines are infinite
            2. get the points from the lines that is a distance (the bend radius) away from the line
            3. choose the point for each line that is in the same direction as the other line
            4. adjust both lines towards each other to go through the closer points
            5. find the intersection of the adjusted lines
        """
        debug = True
        a0, a1 = self.segments[0].xy0, self.segments[0].xy1
        b0, b1 = self.segments[1].xy0, self.segments[1].xy1

        # check that the middle point is shared
        if a1 != b0:
            raise ValueError("Error in TraceCorner.get_arc_seg_intersections(): " + f"the midpoint between the segments ({a0},{a1}) and ({b0},{b1}) isn't shared! " + f"{a1} != {b0}")

        # don't search for the arc segment intersection of two segments with the same(ish) slope
        angle_a = geo.line_segment_angle(a0, a1)
        angle_b = geo.line_segment_angle(b0, b1)
        angle_diff = geo.normalize_angle(angle_b - angle_a)
        if debug:
            print(f"{angle_diff=}")
        if angle_diff < 1/1e6 or angle_diff > (2*np.pi - 1/1e6):
            return None

        # 2. get the bend radius distance points away from each line
        def distance_along_tangent(line: tuple[tuple[float, float], tuple[float, float]], from_point: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
            tangent = geo.get_tangent_line(line, from_point)
            if debug:
                line_slope_intercept = geo.line_two_points_to_slope_intercept(line)
                tangent_slope_intercept = geo.line_two_points_to_slope_intercept(tangent)
                print(f"{line_slope_intercept=}\n{tangent_slope_intercept=}")
            ret0 = geo.distance_along_line((tangent[0], tangent[1]), self.bend_radius, from_point)
            ret1 = geo.distance_along_line((tangent[1], tangent[0]), self.bend_radius, from_point)
            return ret0, ret1
        pnts_dist_adj_a = distance_along_tangent((a0, a1), a1)
        if debug:
            print(f"{pnts_dist_adj_a=}")
        pnts_dist_adj_b = distance_along_tangent((b0, b1), a1)
        if debug:
            print(f"{pnts_dist_adj_b=}")

        # 3. choose the point for each line that is in the same direction as the other line
        def get_point_on_same_side(line: tuple[tuple[float, float], tuple[float, float]],
                                   ref_point: tuple[float, float],
                                   two_points: tuple[tuple[float, float], tuple[float, float]]) -> tuple[float, float]:
            if geo.is_point_on_right(line, two_points[0]):
                if geo.is_point_on_right(line, ref_point):
                    return two_points[0]
                else:
                    return two_points[1]
            else: # point 0 is on the left
                if not geo.is_point_on_right(line, ref_point):
                    return two_points[0]
                else:
                    return two_points[1]
        pnt_dist_adj_a = get_point_on_same_side((a0, a1), b1, pnts_dist_adj_a)
        pnt_dist_adj_b = get_point_on_same_side((b1, b0), a0, pnts_dist_adj_b)
        if debug:
            print(f"{pnt_dist_adj_a=}\n{pnt_dist_adj_b=}")

        # 4. adjust both lines towards each other to go through the closer points
        def get_line_through_point(parallel_line: tuple[tuple[float, float], tuple[float, float]], point: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
            (x1, y1), (x2, y2) = parallel_line
            slope, _ = geo.line_two_points_to_slope_intercept(parallel_line)

            # check for infinity or zero
            if abs(slope) == np.inf:
                if y2 > y1:
                    return point, (point[0], point[1] + 10)
                else:
                    return point, (point[0], point[1] - 10)
            elif slope == 0:
                if x2 > x1:
                    return point, (point[0] + 10, point[1])
                else:
                    return point, (point[0] - 10, point[1])
            
            # get the new line
            return point, (point[0]+10, slope*10 + point[1])
        adj_line_a = get_line_through_point((a0, a1), pnt_dist_adj_a)
        adj_line_b = get_line_through_point((b0, b1), pnt_dist_adj_b)
        if debug:
            print(f"{adj_line_a=}\n{adj_line_b=}")

        # 5. find the intersection of the adjusted lines
        adj_lines_intersection = geo.lines_intersection(adj_line_a, adj_line_b)
        if debug:
            print(f"{adj_lines_intersection=}")

        return adj_lines_intersection

    def get_arc_properties(self) -> tuple[float, float, float, tuple[float, float]]:
        debug = False

        # Get some segment properties
        angle_a = geo.line_segment_angle(self.segments[0].xy0, self.segments[0].xy1)
        angle_b = geo.line_segment_angle(self.segments[1].xy0, self.segments[1].xy1)
        if debug:
            print(f"{angle_a=}\n{angle_b=}")

        # Compute the arc angles.
        # This is more easily done by adjusting the angles such that angle_a = 0.
        angle_b = geo.normalize_angle(angle_b - angle_a)
        angle_diff = angle_b if angle_b < np.pi else angle_b-2*np.pi
        mid_angle = geo.normalize_angle(angle_diff/2 + angle_a)
        if debug:
            print(f"{angle_b=}\n{angle_diff=}\n{angle_diff/2+angle_a=}\n{mid_angle=}")

        # Compute the length of the arc necessary based on
        # the given bend radius and the two angles a and b.
        circle_proportion = abs(angle_diff) / (2*np.pi)
        arc_length = circle_proportion * 2*np.pi * self.bend_radius

        # Get the center for where to place the arc
        arc_center = self.get_arc_seg_intersections()

        return angle_diff, mid_angle, arc_length, arc_center

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


if __name__ == "__main__":
    from SingleTrace import SingleTrace
    from PipeShape import PipeBasicBox
    from units import *

    xy_points = [
        [0, 0],
        [5, 5],
        [10, 5]
    ]
    segments = [
        [0, 1],
        [1, 2]
    ]
    shape = PipeBasicBox(awg2mm(26))
    trace = SingleTrace(xy_points, segments, shape)
    corner = TraceCorner(trace, trace.segments, 1)

    angle_diff, mid_angle, arc_length, arc_center = corner.get_arc_properties()
    print(f"{angle_diff=}\n{mid_angle=}\n{arc_length=}\n{arc_center=}")