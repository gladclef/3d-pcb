from functools import cache
import numpy as np
import vtk

from Geometry.Line import Line
from Geometry.Path import Path
import Geometry.geometry_tools as geo
from Geometry.LineSegment import LineSegment

class PathCorner:
    def __init__(self, parent: Path, segments: tuple[LineSegment, LineSegment], bend_radius: float):
        """
        Parameters
        ----------
        parent : Path
            The parent instance that contains this instance.
        segments : tuple[LineSegment, LineSegment]
            The two segments that this corner is found on. They should
            share one of the end points.
        bend_radius : float
            How large of a radius the bend should be.
        """
        # verify that the segments share an end point, and that the shared point is in the middle
        if segments[0].xy2 == segments[1].xy1:
            if segments[0].xy1 == segments[1].xy2:
                raise ValueError("The segments must share only one common end point.")
        elif segments[0].xy1 == segments[1].xy2:
            segments = [segments[1], segments[0]]
        else:
            raise ValueError("The segments must share a common end point.")
        
        # check for a valid bend radius
        if bend_radius < parent.shape.radius:
            raise RuntimeError("Error in PathCorner.min_corner_increment_angle: " + f"bend_radius {bend_radius} is less than the pipe shape radius {parent.shape.radius}!")

        self.parent = parent
        self.segments = segments
        self._bend_radius = bend_radius
        self.debug = False

    @property
    def bend_radius(self):
        return self._bend_radius
    
    @bend_radius.setter
    def bend_radius(self, val: float):
        self._bend_radius = val

    def get_arc_seg_intersections(self):
        """
        Find the location where the arc intersects both line segments.

        To do this we:
            1. assume both lines are infinitely long
            2. get the points from the lines that is a distance (the bend radius) away from the line
            3. choose the point for each line that is in the same direction as the other line
            4. adjust both lines towards each other to go through the closer points
            5. find the intersection of the adjusted lines
        """
        seg1 = self.segments[0]
        seg2 = self.segments[1]
        line1 = Line.from_two_points(seg1.xy1, seg1.xy2)
        line2 = Line.from_two_points(seg2.xy1, seg2.xy2)

        # check that the middle point is shared
        if seg1.xy2 != seg2.xy1:
            raise ValueError("Error in PathCorner.get_arc_seg_intersections(): " + f"the midpoint between the segments ({seg1.xy2}) and ({seg2.xy1}) isn't shared!")

        # don't search for the arc segment intersection of two segments with the same(ish) slope
        angle_a = line1.angle
        angle_b = line2.angle
        angle_diff = geo.normalize_angle(angle_b - angle_a)
        if self.debug:
            print(f"{angle_diff=}")
        if angle_diff < geo.ZERO_THRESH or angle_diff > (2*np.pi - geo.ZERO_THRESH):
            tangent = line1.get_tangent_line(seg1.xy2)
            return tangent.distance_along_line(self.bend_radius, seg1.xy2)

        # 2. get the bend radius distance points away from each line
        def distance_along_tangent(line: Line, from_point: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
            tangent = line.get_tangent_line(from_point)
            if self.debug:
                line_slope_intercept = line.slope, line.y_intercept
                tangent_slope_intercept = tangent.slope, tangent.y_intercept
                print(f"{line_slope_intercept=}\n{tangent_slope_intercept=}")
            ret0 = tangent.distance_along_line(self.bend_radius, from_point)
            ret1 = tangent.reversed().distance_along_line(self.bend_radius, from_point)
            return ret0, ret1
        pnts_dist_adj_a = distance_along_tangent(line1, seg1.xy2)
        if self.debug:
            print(f"{pnts_dist_adj_a=}")
        pnts_dist_adj_b = distance_along_tangent(line2, seg1.xy2)
        if self.debug:
            print(f"{pnts_dist_adj_b=}")

        # 3. choose the point for each line that is in the same direction as the other line
        def get_point_on_same_side(line: Line,
                                   ref_point: tuple[float, float],
                                   two_points: tuple[tuple[float, float], tuple[float, float]]) -> tuple[float, float]:
            if line.is_point_on_right(two_points[0]):
                if line.is_point_on_right(ref_point):
                    return two_points[0]
                else:
                    return two_points[1]
            else: # point 0 is on the left
                if not line.is_point_on_right(ref_point):
                    return two_points[0]
                else:
                    return two_points[1]
        pnt_dist_adj_a = get_point_on_same_side(line1, seg2.xy2, pnts_dist_adj_a)
        pnt_dist_adj_b = get_point_on_same_side(line2, seg1.xy1, pnts_dist_adj_b)
        if self.debug:
            print(f"{pnt_dist_adj_a=}\n{pnt_dist_adj_b=}")

        # 4. adjust both lines towards each other to go through the closer points
        adj_line_a = Line.from_angle_point(line1.angle, pnt_dist_adj_a)
        adj_line_b = Line.from_angle_point(line2.angle, pnt_dist_adj_b)
        if self.debug:
            print(f"{adj_line_a=}\n{adj_line_b=}")

        # 5. find the intersection of the adjusted lines
        adj_lines_intersection = adj_line_a.intersection(adj_line_b)
        if self.debug:
            print(f"{adj_lines_intersection=}")

        return adj_lines_intersection

    def get_arc_properties(self) -> tuple[float, float, float, tuple[float, float]]:
        """ Get common properties about the arc that this path corner follows.

        Returns
        -------
        angle_diff: float
            The difference between the angle of the first and second segments.
            This value plus the angle of the first segment will equal the angle
            of second segment (angle_seg1 + angle_diff == angle_seg2).
        mid_angle: float
            The angle that is halfway between the angle of the first and second segments.
        arc_length: float
            The total length of the arc. Note that this is the length for a perfect
            arc. The actual length will be different due to the path corner being
            composed of multiple increments.
        arc_center: tuple[float, float]
            The x,y point at the center of the circle that the arc is a part of.
        """

        # Get some segment properties
        angle_a = self.segments[0].angle
        angle_b = self.segments[1].angle
        if self.debug:
            print(f"{angle_a=}\n{angle_b=}")

        # Compute the arc angles.
        # This is more easily done by adjusting the angles such that angle_a = 0.
        angle_b = geo.normalize_angle(angle_b - angle_a)
        angle_diff = angle_b if angle_b < np.pi else angle_b-2*np.pi
        mid_angle = geo.normalize_angle(angle_diff/2 + angle_a)
        if self.debug:
            print(f"{angle_b=}\n{angle_diff=}\n{angle_diff/2+angle_a=}\n{mid_angle=}")

        # Compute the length of the arc necessary based on
        # the given bend radius and the two angles a and b.
        circle_proportion = abs(angle_diff) / (2*np.pi)
        arc_length = circle_proportion * 2*np.pi * self.bend_radius

        # Get the center for where to place the arc
        arc_center = self.get_arc_seg_intersections()

        return angle_diff, mid_angle, arc_length, arc_center

    @cache
    def get_center_line_points(self):
        """
        Get the points that define the path for the path corner.

        These are the points along the "center line" of the path corner,
        as pairs of (xy, θ), where:
            - xy is the x,y location of the point
            - θ is the direction that point is facing (0-2π)
        
        The first point in the path corner corresponds to the end of the
        first segment, and the last point corresponds to the start of the
        second segment. There will therefore always be at least two points.

        Returns
        -------
        xypoint_angle_pairs: list[tuple[tuple[float, float], float]]
            Pairs of (xy, θ), one per point on the center line.
        """
        corner_increment_angle = np.deg2rad(5)
        angle_diff, mid_angle, arc_length, arc_center = self.get_arc_properties()

        # Determine the number of increments to be used.
        # Note that the start of the first increment and the end of the last
        # increment is where the segment ends will be.
        nincrements = int(abs(angle_diff) / (corner_increment_angle*2 - np.deg2rad(1))) - 1
        nincrements = max(nincrements, 1)

        # Find the tangent from the arc center to the path corner (aka the end of the first segment).
        seg_tangent = self.segments[0].get_tangent_line(self.segments[0].xy2)
        tan_tst_pnt = seg_tangent.distance_along_line(1.0, self.segments[0].xy2)
        if self.segments[0].is_point_on_right(arc_center):
            if self.segments[0].is_point_on_right(tan_tst_pnt):
                seg_tangent = seg_tangent.reversed()
                tan_tst_pnt2 = seg_tangent.distance_along_line(1.0, self.segments[0].xy2)
                assert not self.segments[0].is_point_on_right(tan_tst_pnt2)
        else:
            if not self.segments[0].is_point_on_right(tan_tst_pnt):
                seg_tangent = seg_tangent.reversed()
                tan_tst_pnt2 = seg_tangent.distance_along_line(1.0, self.segments[0].xy2)
                assert self.segments[0].is_point_on_right(tan_tst_pnt2)

        # Get the angle for the start of the first increment.
        first_tan_angle = seg_tangent.angle

        # Build out the list of increments.
        first_angle = self.segments[0].angle
        inc_angle = angle_diff / nincrements
        ret: list[tuple[tuple[float, float], float]] = []
        for i in range(nincrements+1):
            tan_angle = first_tan_angle + i*inc_angle
            pnt_xy = np.cos(tan_angle)*self.bend_radius, np.sin(tan_angle)*self.bend_radius
            pnt_xy = pnt_xy[0] + arc_center[0], pnt_xy[1] + arc_center[1]
            pnt_angle = first_angle + i*inc_angle
            ret.append((pnt_xy, pnt_angle))
        
        return ret

    def draw_corner(self):
        """ For debugging: draw the corner center line to be followed during shape generation. """
        import matplotlib.pyplot as plt
        
        angle_diff, mid_angle, arc_length, arc_center = self.get_arc_properties()
        center_line_points = self.get_center_line_points()

        # create the plot
        fig, ax = plt.subplots(figsize=(10,10))

        # draw the segments
        dx0, dy0 = self.segments[0].x2 - self.segments[0].x1, self.segments[0].y2 - self.segments[0].y1
        ax.arrow(*self.segments[0].xy1, dx0, dy0, color="tab:blue")
        dx1, dy1 = self.segments[1].x2 - self.segments[1].x1, self.segments[1].y2 - self.segments[1].y1
        ax.arrow(*self.segments[1].xy1, dx1, dy1, color="tab:blue")

        # draw the center point
        ax.add_patch(plt.Circle(arc_center, self.bend_radius/10, color="tab:purple"))

        # draw the arc
        for pnt_xy, pnt_ang in center_line_points:
            ax.arrow(*pnt_xy, np.cos(pnt_ang)*self.bend_radius/2, np.sin(pnt_ang)*self.bend_radius/2, color="black")
        for i in range(len(center_line_points)-1):
            pnt0_xy, pnt0_ang = center_line_points[i]
            pnt1_xy, pnt1_ang = center_line_points[i+1]
            ax.plot([pnt0_xy[0], pnt1_xy[0]], [pnt0_xy[1], pnt1_xy[1]], color="tab:orange")
        
        # show the plot
        ax.set_aspect('equal')
        ax.set_xlim(-3.5, 3.5)
        ax.set_ylim(-3.5, 3.5)
        plt.show(block=False)


if __name__ == "__main__":
    import json

    import matplotlib.pyplot as plt

    from Trace.SingleTrace import SingleTrace
    from Trace.PipeShape import DEFAULT_PIPE_SHAPE
    from tool.units import *

    vidframes_out = open("output/renders/input.txt", "w")
    unittests_out = None #open("test/test_PathCorners_regression_vals.json", "w")

    if unittests_out is not None:
        unittests_out.write("[\n")

    for i, deg100 in enumerate(range(14500, -14500, -25)):
        deg = deg100 / 100
        rad = np.deg2rad(deg)

        xy_points = [
            [-3.5, 0],
            [0, 0],
            [0+np.cos(rad)*3.5, 0+np.sin(rad)*3.5]
        ]
        segments = [
            [0, 1],
            [1, 2]
        ]
        shape = DEFAULT_PIPE_SHAPE(awg2mm(26))
        trace = SingleTrace(xy_points, segments, shape)
        corner = PathCorner(trace, trace.segments, 1)

        if vidframes_out is not None:
            corner.draw_corner()
            plt.savefig(f"output/renders/{i}_{deg}.png")
            plt.close('all')
            vidframes_out.write(f"file '{i}_{deg}.png'\n")
            vidframes_out.write(f"duration {1/30}\n")
            # ffmpeg command: something like "ffmpeg -f concat -i 'input.txt' -c:v -r 20 out.mp4"
        
        if unittests_out is not None:
            angle_diff, mid_angle, arc_length, arc_center = corner.get_arc_properties()
            center_line_points = corner.get_center_line_points()
            regression_val = {
                "deg100": deg100,
                "deg": deg,
                "rad": rad,
                "angle_diff": angle_diff,
                "mid_angle": mid_angle,
                "arc_length": arc_length,
                "arc_center": arc_center,
                "center_line_points": center_line_points
            }
            regression_str = "" if i == 0 else ",\n"
            regression_str += json.dumps(regression_val).replace("\n", " ")
            unittests_out.write(regression_str)
    
    if vidframes_out is not None:
        vidframes_out.close()
    if unittests_out is not None:
        unittests_out.write("\n]")
        unittests_out.close()